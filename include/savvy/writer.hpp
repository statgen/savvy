/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_BCF_WRITER_HPP
#define LIBSAVVY_BCF_WRITER_HPP

#include <unistd.h>

#include "file.hpp"
#include "typed_value.hpp"
#include "utility.hpp"
#include "compressed_vector.hpp"
#include "region.hpp"
#include "s1r.hpp"
#include "pbwt.hpp"


#include <shrinkwrap/zstd.hpp>
#include <shrinkwrap/gz.hpp>

#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <chrono>
#include <list>
#include <cstdint>
#include <type_traits>

namespace savvy
{
  namespace v2
  {
//    template <typename DestT, typename InT, typename OutT>
//    void compress_sparse_offsets(InT in, InT in_end, OutT out, std::size_t stride)
//    {
//      std::size_t last_off = 0;
//      //char* char_p = (char*)out;
//      for (auto it = in; it != in_end; ++it)
//      {
//        DestT off((*it) - last_off);
//        *(out++) = off;
//        //*out = off;
//        //char_p += stride;
//        //out = (OutT)char_p;
//        last_off = (*it) + 1;
//      }
//    }

    class writer : public file
    {
    public:
      static const int default_compression_level = 3;
      static const int default_block_size = 2048;
    private:
      format file_format_;
      std::mt19937_64 rng_;
      std::unique_ptr<std::streambuf> output_buf_;
      std::string file_path_;
      std::ostream ofs_;
      std::array<std::uint8_t, 16> uuid_;
      std::size_t n_samples_ = 0;
      std::vector<char> serialized_buf_;
      phasing phasing_ = phasing::unknown;

      // Data members to support indexing
      std::unique_ptr<s1r::writer> index_file_;
      std::string current_chromosome_;
      std::size_t block_size_ = default_block_size;
      std::size_t record_count_ = 0;
      std::size_t record_count_in_block_ = 0;
      std::uint32_t current_block_min_ = std::numeric_limits<std::uint32_t>::max();
      std::uint32_t current_block_max_ = 0;

      ::savvy::internal::pbwt_sort_context sort_context_;
    private:
      static std::filebuf *create_std_filebuf(const std::string& file_path, std::ios::openmode mode);

      static std::unique_ptr<std::streambuf> create_out_streambuf(const std::string& file_path, format file_format, std::uint8_t compression_level);

    public:
      writer(const std::string& file_path, file::format file_format, std::vector<std::pair<std::string, std::string>> headers, const std::vector<std::string>& ids, std::uint8_t compression_level = default_compression_level, bool create_index = true);

      ~writer();

      void set_block_size(std::uint16_t bs);

      bool good() const { return ofs_.good(); }
      operator bool() const { return ofs_.good(); }
      //bool bad() const { return ofs_.bad(); }

      writer& write(const variant& r);
      writer& operator<<(const variant& v) { return write(v); }

    private:
      writer& write_vcf(const variant& r);
      void write_header(std::vector<std::pair<std::string, std::string>>& headers, const std::vector<std::string>& ids);

      bool serialize_vcf_shared(const site_info& s);
      bool serialize_vcf_indiv(const variant& v, phasing phased);
    };


    //================================================================//
    // Writer definitions
    inline
    std::filebuf* writer::create_std_filebuf(const std::string& file_path, std::ios::openmode mode)
    {
      std::filebuf* ret = new std::filebuf();
      ret->open(file_path.c_str(), mode);
      return ret;
    }

    inline
    std::unique_ptr<std::streambuf> writer::create_out_streambuf(const std::string& file_path, format file_fmt, std::uint8_t compression_level)
    {
      if (compression_level > 0)
      {
        if (file_fmt == format::sav2 || file_fmt == format::sav1)
          return std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path, compression_level));
        else
          return std::unique_ptr<std::streambuf>(new shrinkwrap::bgzf::obuf(file_path));  //, compression_level)); TODO: Add compression level
      }
      else
      {
        return std::unique_ptr<std::streambuf>(create_std_filebuf(file_path, std::ios::binary | std::ios::out));
      }
    }

    inline
    writer::writer(const std::string& file_path, file::format file_format, std::vector<std::pair<std::string, std::string>> headers, const std::vector<std::string>& ids, std::uint8_t compression_level, bool create_index) :
      file_format_(file_format),
      rng_(std::chrono::high_resolution_clock::now().time_since_epoch().count() ^ std::clock() ^ (std::uint64_t) this),
      output_buf_(create_out_streambuf(file_path, file_format_, compression_level)), //opts.compression == compression_type::zstd ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path)) : std::unique_ptr<std::streambuf>(new std::filebuf(file_path, std::ios::binary))),
      ofs_(output_buf_.get()),
      //samples_(samples_beg, samples_end),
      file_path_(file_path),
      uuid_(::savvy::detail::gen_uuid(rng_))
    {
      if (file_format_ == format::sav1)
      {
        fprintf(stderr, "Error: Writing SAV v1 format not supported.\n");
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
        return;
      }

      // TODO: Use mkstemp when shrinkwrap supports FILE*
      if (file_format_ == format::sav2 && create_index) // TODO: Check if zstd is enabled.
      {
        std::string idx_path = "/tmp/tmpfileXXXXXX";
        int tmp_fd = mkstemp(&idx_path[0]);
        if (!tmp_fd)
        {
          std::cerr << "Error: could not open temp file for s1r index (" << idx_path << ")" << std::endl;
        }
        else
        {
          index_file_ = ::savvy::detail::make_unique<s1r::writer>(idx_path, uuid_);
          std::remove(idx_path.c_str());
          ::close(tmp_fd);
        }
      }

      write_header(headers, ids);
    }

    inline
    writer::~writer()
    {
      // TODO: This is only a temp solution.
      if (index_file_)
      {
        if (record_count_in_block_)
        {
          auto file_pos = std::uint64_t(ofs_.tellp());
          if (record_count_in_block_ > 0x10000) // Max records per block: 64*1024
          {
            assert(!"Too many records in zstd frame to be indexed!");
          }

          if (file_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
          {
            assert(!"File size to large to be indexed!");
          }

          s1r::entry e(current_block_min_, current_block_max_, (file_pos << 16) | std::uint16_t(record_count_in_block_ - 1));
          index_file_->write(current_chromosome_, e);
        }

        ofs_.flush();
        auto idx_fs = index_file_->close();
        std::fstream ofs(file_path_, std::ios::out | std::ios::binary | std::ios::app); // TODO: THIS SEEMS DANGEROUS. Store FILE* when creating zstd stream and use instead of opening new descriptor.
        std::int64_t p = ofs.tellp();
        if (!::savvy::detail::append_skippable_zstd_frame(idx_fs, ofs))
        {
          ofs_.setstate(ofs_.rdstate() | std::ios::badbit); // TODO: Use linkat or send file (see https://stackoverflow.com/a/25154505/1034772)
          std::cerr << "Error: index file too big for skippable zstd frame" << std::endl;
        }
      }
    }

    inline
    void writer::set_block_size(std::uint16_t bs)
    {
      block_size_ = bs;
    }

    inline
    writer& writer::write_vcf(const variant& r)
    {
      if (!serialize_vcf_shared(r) || !serialize_vcf_indiv(r, phasing_))
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit);

      return *this;
    }

    inline
    writer& writer::write(const variant& r)
    {
      if (file_format_ == format::vcf)
        return write_vcf(r);

      bool is_bcf = file_format_ == format::bcf; // TODO: ...
      bool flushed = false;

      if (block_size_ != 0 && (block_size_ <= record_count_in_block_ || r.chrom() != current_chromosome_)) // TODO: this needs to be fixed to support variable block size
      {
        if (index_file_ && record_count_in_block_)
        {
          auto file_pos = std::uint64_t(ofs_.tellp());
          if (record_count_in_block_ > 0x10000) // Max records per block: 64*1024
          {
            assert(!"Too many records in zstd frame to be indexed!");
            ofs_.setstate(std::ios::badbit);
          }

          if (file_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
          {
            assert(!"File size to large to be indexed!");
            ofs_.setstate(std::ios::badbit);
          }

          s1r::entry e(current_block_min_, current_block_max_, (file_pos << 16) | std::uint16_t(record_count_in_block_ - 1));
          index_file_->write(current_chromosome_, e);
        }
        ofs_.flush();
        current_chromosome_ = r.chrom();
        record_count_in_block_ = 0;
        current_block_min_ = std::numeric_limits<std::uint32_t>::max();
        current_block_max_ = 0;

        sort_context_.reset();
        flushed = true;
      }

//        std::vector<std::string> pbwt_info_flags;
//        pbwt_info_flags.reserve(fmt_to_pbwt_context_.size());
//        for (auto it = r.format_fields().begin(); it != r.format_fields().end(); ++it)
//        {
//          if (!it->second.is_sparse())
//          {
//            auto res = fmt_to_pbwt_context_.find(it->first);
//            for (auto jt = res; jt != r.format_fields().end(); ++jt)
//            {
//              if (jt->second.get().ploidy == it->second.size() / n_samples_)
//              {
//                pbwt_info_flags.emplace_back(jt->second.get().id);
//                break;
//              }
//            }
//          }
//        }

      std::uint32_t shared_sz, indiv_sz;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // Determine which fields sort
      std::list<std::pair<std::string, typed_value>> extra_info_fields;
      if (flushed && sort_context_.format_contexts.size())
      {
        extra_info_fields.emplace_back("_PBWT_RESET", typed_value(std::int8_t(1)));
      }

      std::size_t n_fmt = 0;
      std::vector<::savvy::internal::pbwt_sort_format_context*> pbwt_format_pointers;
      pbwt_format_pointers.reserve(r.format_fields().size());
      for (auto it = r.format_fields().begin(); it != r.format_fields().end(); ++it)
      {
        pbwt_format_pointers.emplace_back(nullptr);
        if (!it->second.is_sparse())
        {
          auto f = sort_context_.field_to_format_contexts.find(it->first);
          if (f != sort_context_.field_to_format_contexts.end())
          {
            if (f->second->sort_map.size() == it->second.size() || f->second->sort_map.empty())
            {
              pbwt_format_pointers.back() = f->second;
              extra_info_fields.emplace_back(f->second->id, typed_value(std::int8_t(1)));
            }
          }
        }
        if (file_format_ == format::sav2 || it->first != "PH")
          ++n_fmt;
      }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

      serialized_buf_.clear();
      serialized_buf_.reserve(24);

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // Serialize shared data
      if (!site_info::serialize(r, std::back_inserter(serialized_buf_), dict_, n_samples_, n_fmt, extra_info_fields, sort_context_.format_contexts.size() ? flushed : false))
      {
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
        return *this;
      }

      if (serialized_buf_.size() > std::numeric_limits<std::uint32_t>::max())
      {
        fprintf(stderr, "Error: shared data too big\n");
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit); // TODO: Maybed some of these should be fail instead of bad.
        return *this;
      }

      shared_sz = serialized_buf_.size();
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // Serialize individual data
      if (!variant::serialize(r, std::back_inserter(serialized_buf_),
        dict_, n_samples_, is_bcf, phasing_,
        sort_context_, pbwt_format_pointers))
      {
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
        return *this;
      }


      if (serialized_buf_.size() - shared_sz > std::numeric_limits<std::uint32_t>::max())
      {
        std::fprintf(stderr, "Error: individual data too big\n");
        ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
        return *this;
      }

      indiv_sz = serialized_buf_.size() - shared_sz;
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


      ofs_.write((char *) &shared_sz, sizeof(shared_sz));
      ofs_.write((char *) &indiv_sz, sizeof(indiv_sz));
      ofs_.write(serialized_buf_.data(), serialized_buf_.size());

      current_block_min_ = std::min(current_block_min_, std::uint32_t(r.pos()));
      std::size_t max_alt_size = 0;
      for (auto it = r.alts().begin(); it != r.alts().end(); ++it)
        max_alt_size = std::max(max_alt_size, it->size());
      current_block_max_ = std::max(current_block_max_, std::uint32_t(r.pos() + std::max(r.ref().size(), max_alt_size)) - 1);
      ++record_count_in_block_;
      ++record_count_;


      return *this;
    }

    inline
    void writer::write_header(std::vector<std::pair<std::string, std::string>>& headers, const std::vector<std::string>& ids)
    {
      std::string magic = {'S', 'A', 'V', '\x02', '\x00'};

      dict_.str_to_int[dictionary::id]["PASS"] = dict_.entries[dictionary::id].size();
      dict_.entries[dictionary::id].emplace_back(dictionary::entry{"PASS", "", 0});

      bool gt_present{}, ph_present{};

      std::uint32_t header_block_sz = 0;
      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        auto hval = parse_header_value(it->second);
        if (!hval.id.empty() && hval.id != "PASS")
        {
          int which_dict = -1;
          if (it->first == "contig") which_dict = dictionary::contig;
          else if (it->first == "INFO" || it->first == "FILTER" || it->first == "FORMAT") which_dict = dictionary::id;
          else if (it->first == "SAMPLE") which_dict = dictionary::sample;

          if (which_dict >= 0 && dict_.str_to_int[which_dict].find(hval.id) == dict_.str_to_int[which_dict].end())
          {
            dictionary::entry e;
            e.id = hval.id;
            e.number = std::atoi(hval.number.c_str()); // TODO: handle special character values.
            if (hval.type == "Integer")
              e.type = typed_value::int32;
            else if (hval.type == "Float")
              e.type = typed_value::real;
            else if (hval.type == "String")
              e.type = typed_value::str;

            if (!hval.idx.empty())
            {
              std::size_t idx = std::atoi(hval.idx.c_str());
              dict_.entries[which_dict].resize(idx, {"DELETED", "", 0});
            }

            dict_.str_to_int[which_dict][hval.id] = dict_.entries[which_dict].size();
            dict_.entries[which_dict].emplace_back(std::move(e));
          }
        }

        if (it->first == "INFO")
        {
          if (hval.id.substr(0, 10) == "_PBWT_SORT")
          {
            ::savvy::internal::pbwt_sort_format_context ctx;
            ctx.format = parse_header_sub_field(it->second, "Format");
            ctx.id = hval.id;

            auto insert_it = sort_context_.format_contexts.insert(std::make_pair(std::string(hval.id), std::move(ctx)));
            sort_context_.field_to_format_contexts.insert(std::make_pair(insert_it.first->second.format, &(insert_it.first->second)));
          }
          // TODO: set _PBWT_RESET flag.
        }
        else if (it->first == "phasing")
        {
          if (it->second == "none")
            phasing_ = phasing::none;
          else if (it->second == "partial")
            phasing_ = phasing::partial;
          else if (it->second == "phased" || it->second == "full")
            phasing_ = phasing::phased;
        }
        else if (it->first == "FORMAT")
        {
          if (hval.id == "GT")
          {
            if (file_format_ == format::sav2)
            {
              if (hval.type == "String")
                it->second = "<ID=GT, Type=Integer, Number=G, Description=\"Genotype\">";
            }
            else
            {
              if (hval.type != "String")
                it->second = "<ID=GT, Type=String, Number=1, Description=\"Genotype\">";
            }

            gt_present = true;
          }
          else if (hval.id == "PH")
          {
            ph_present = true;
          }
        }

        header_block_sz += it->first.size();
        header_block_sz += it->second.size();
        header_block_sz += 4;
      }

      if ((phasing_ == phasing::unknown || phasing_ == phasing::partial) && gt_present && !ph_present) // TODO: potentially make unkkown = none
      {
        headers.emplace_back("FORMAT", "<ID=PH, Type=Integer, Number=., Description=\"Genotype phase\">");

        header_block_sz += headers.back().first.size();
        header_block_sz += headers.back().second.size();
        header_block_sz += 4;

        dict_.str_to_int[dictionary::id]["PH"] = dict_.entries[dictionary::id].size();
        dict_.entries[dictionary::id].emplace_back(dictionary::entry{"PH", ".", typed_value::int8});
      }

      std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
      header_block_sz += column_names.size();

      for (auto it = ids.begin(); it != ids.end(); ++it)
      {
        header_block_sz += 1 + it->size();
      }

      header_block_sz += 2; //new line and null

      if (file_format_ != format::vcf)
      {
        if (file_format_ == format::bcf)
          magic = {'B', 'C', 'F', '\x02', '\x02'};
        ofs_.write(magic.data(), magic.size());
        ofs_.write((char *) (&header_block_sz), sizeof(header_block_sz));
      }

      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        ofs_.write("##", 2);
        ofs_.write(it->first.data(), it->first.size());
        ofs_.write("=", 1);
        ofs_.write(it->second.data(), it->second.size());
        ofs_.write("\n", 1);
      }

      ofs_.write(column_names.data(), column_names.size());
      for (auto it = ids.begin(); it != ids.end(); ++it)
      {
        ofs_.write("\t", 1);
        ofs_.write(it->data(), it->size());
      }
      n_samples_ = ids.size();

      ofs_.put('\n');
      if (file_format_ != format::vcf)
        ofs_.put('\0');
    }

    inline
    bool writer::serialize_vcf_shared(const site_info& s)
    {
      (ofs_) << s.chrom_
        << "\t" << s.pos_
        << "\t" << std::string(s.id_.size() ? s.id_ : ".")
        << "\t" << s.ref_;

      if (s.alts_.empty())
      {
        ofs_ << "\t.";
      }
      else
      {
        ofs_ << "\t" << s.alts_.front();
        for (auto it = s.alts_.begin() + 1; it != s.alts_.end(); ++it)
          ofs_ << "," << *it;
      }

      if (std::isnan(s.qual_))
        ofs_ << "\t.";
      else
        ofs_ << "\t" << s.qual_;

      if (s.filters_.empty())
      {
        ofs_ << "\t.";
      }
      else
      {
        ofs_ << "\t" << s.filters_.front();
        for (auto it = s.filters_.begin() + 1; it != s.filters_.end(); ++it)
          ofs_ << "," << *it;
      }

      if (s.info_.empty())
      {
        ofs_ << "\t.";
      }
      else
      {
        ofs_ << "\t" << s.info_.front().first << "=" << s.info_.front().second;
        for (auto it = s.info_.begin() + 1; it != s.info_.end(); ++it)
        {
          ofs_ << ";" << it->first << "=" << it->second;
        }
      }

      return ofs_.good();
    }

    inline
    bool writer::serialize_vcf_indiv(const savvy::v2::variant& v, phasing phased)
    {
      std::vector<const typed_value*> typed_value_ptrs(v.format_fields_.size());
      std::vector<std::size_t> strides(v.format_fields_.size());
      std::vector<char> delims(v.format_fields_.size());
      std::list<typed_value> dense_copies;
      std::int8_t* ph_ptr = nullptr;
      for (std::size_t i = 0; i < v.format_fields_.size(); ++i)
      {
        assert(n_samples_); // TODO

        if (v.format_fields_[i].first == "PH")
        {
          assert(i == 1); // TODO: return error
          ph_ptr = (std::int8_t*)v.format_fields_[i].second.val_ptr_;
          continue;
        }

        ofs_ << (i == 0 ? "\t" : ":") << v.format_fields_[i].first;

        strides[i] = (n_samples_ ? v.format_fields_[i].second.size() / n_samples_ : 0);

        if (v.format_fields_[i].first == "GT")
          delims[i] = (phased == phasing::phased ? '|' : '/');
        else
          delims[i] = ',';

        if (v.format_fields_[i].second.is_sparse())
        {
          dense_copies.emplace_back();
          v.format_fields_[i].second.copy_as_dense(dense_copies.back());
          typed_value_ptrs[i] = &dense_copies.back();
        }
        else
        {
          typed_value_ptrs[i] = &v.format_fields_[i].second;
        }
      }


      if (ph_ptr)
      {
        std::size_t ph_stride = strides[0] - 1;
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          ofs_.put('\t');
          for (std::size_t k = 0; k < strides[0]; ++k)
          {
            if (k > 0)
              ofs_.put(ph_ptr[i * ph_stride + k - 1] ? '|' : '/'); // TODO: allow for PH
            typed_value_ptrs[0]->serialize_vcf(i * strides[0] + k, ofs_);
          }

          for (std::size_t j = 2; j < v.format_fields_.size(); ++j)
          {
            ofs_.put(':');
            for (std::size_t k = 0; k < strides[j]; ++k)
            {
              if (k > 0)
                ofs_.put(delims[j]);
              typed_value_ptrs[j]->serialize_vcf(i * strides[j] + k, ofs_);
            }
          }
        }
      }
      else
      {
        for (std::size_t i = 0; i < n_samples_; ++i)
        {
          for (std::size_t j = 0; j < v.format_fields_.size(); ++j)
          {
            ofs_.put(j > 0 ? ':' : '\t');
            for (std::size_t k = 0; k < strides[j]; ++k)
            {
              if (k > 0)
                ofs_.put(delims[j]); // TODO: allow for PH
              typed_value_ptrs[j]->serialize_vcf(i * strides[j] + k, ofs_);
            }
          }
        }
      }

      ofs_.put('\n');

      return ofs_.good();
    }
    //================================================================//




//    template <typename T>
//    void variant::serialize_format(std::vector<char>& dest, const std::vector<T>& geno, bool sparse)
//    {
//
//    }
//
//    template <typename T>
//    void variant::serialize_format(std::vector<char>& dest, const compressed_vector<T>& geno)
//    {
//      std::size_t offset_max = 0;
//      std::size_t last_off = 0;
//      for (auto it = geno.begin(); it != geno.end(); ++it)
//      {
//        std::size_t off = it.offset() - last_off;
//        last_off = it.offset() + 1;
//        if (off > offset_max)
//          offset_max = off;
//      }
//
//      dest.clear();
//      dest.reserve(1 + 8);
//
//
//      //std::uint8_t off_type = bcf::int_type(geno.size() ? (std::int64_t)geno.size() - 1 : 0);
//      std::uint8_t off_type = bcf::int_type((std::int64_t)offset_max);
//      std::uint8_t val_type = sav2::typed_value_old::type_code<T>();
//
//      int off_size = 1u << bcf_type_shift[off_type];
//      int val_size = 1u << bcf_type_shift[val_type];
//      int pair_size = off_size + val_size;
//
//      if (val_size * geno.size() < pair_size * geno.non_zero_size())
//      {
//        std::uint8_t type_byte = val_size;
//        if (geno.size() >= 15u)
//          type_byte = (15u << 4u) | type_byte;
//        else
//          type_byte = (geno.size() << 4u) | type_byte;
//        dest.emplace_back(type_byte);
//
//        if (geno.size() >= 15u)
//        {
//          bcf::serialize_typed_scalar(std::back_inserter(dest), (std::int32_t)geno.size());
//        }
//
//        std::size_t data_pos = dest.size();
//        dest.resize(data_pos + (geno.size() * val_size));
//        if (val_type == 1)
//        {
//          auto out = (std::int8_t*)(dest.data() + data_pos);
//          for (auto it = geno.begin(); it != geno.end(); ++it)
//          {
//            out[it.offset()] = *it;
//          }
//        }
//        else if (val_type == 2)
//        {
//          auto out = (std::int16_t*)(dest.data() + data_pos);
//          for (auto it = geno.begin(); it != geno.end(); ++it)
//          {
//            out[it.offset()] = *it;
//          }
//        }
//        else if (val_type == 3)
//        {
//          auto out = (std::int32_t*)(dest.data() + data_pos);
//          for (auto it = geno.begin(); it != geno.end(); ++it)
//          {
//            out[it.offset()] = *it;
//          }
//        }
//        else if (val_type == 4)
//        {
//          auto out = (std::int64_t*)(dest.data() + data_pos);
//          for (auto it = geno.begin(); it != geno.end(); ++it)
//          {
//            out[it.offset()] = *it;
//          }
//        }
//        else if (val_type == 5)
//        {
//          auto out = (float*)(dest.data() + data_pos);
//          for (auto it = geno.begin(); it != geno.end(); ++it)
//          {
//            out[it.offset()] = *it;
//          }
//        }
////          else if (val_type == 6)
////          {
////            auto out = (double*)(dest.data() + data_pos);
////            for (auto it = geno.begin(); it != geno.end(); ++it)
////            {
////              out[it.offset()] = *it;
////            }
////          }
//        else
//          throw std::runtime_error("Unsupported FORMAT type!");
//      }
//      else
//      {
//        std::uint8_t type_byte = typed_value_old::sparse;
//        if (geno.size() >= 15u)
//          type_byte = (15u << 4u) | type_byte;
//        else
//          type_byte = (geno.size() << 4u) | type_byte;
//        dest.emplace_back(type_byte);
//
//        if (geno.size() >= 15u)
//        {
//          bcf::serialize_typed_scalar(std::back_inserter(dest), (std::int32_t)geno.size());
//        }
//
//        dest.emplace_back((off_type << 4u) | val_type);
//        bcf::serialize_typed_scalar(std::back_inserter(dest), (std::int32_t)geno.non_zero_size());
//
//        std::size_t data_pos = dest.size();
//        dest.resize(data_pos + (geno.non_zero_size() * pair_size));
//
//
//        if (off_type == 1)
//          compress_sparse_offsets<std::int8_t>(geno.index_data(), geno.index_data() + geno.non_zero_size(), (std::int8_t*) (dest.data() + data_pos), pair_size);
//        else if (off_type == 2)
//          compress_sparse_offsets<std::int16_t>(geno.index_data(), geno.index_data() + geno.non_zero_size(), (std::int16_t*) (dest.data() + data_pos), pair_size);
//        else if (off_type == 3)
//          compress_sparse_offsets<std::int32_t>(geno.index_data(), geno.index_data() + geno.non_zero_size(), (std::int32_t*) (dest.data() + data_pos), pair_size);
//        else if (off_type == 4)
//          compress_sparse_offsets<std::int64_t>(geno.index_data(), geno.index_data() + geno.non_zero_size(), (std::int64_t*) (dest.data() + data_pos), pair_size);
//
//        data_pos += off_size * geno.non_zero_size();
//
//        if (val_type == 1)
//          std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (std::int8_t*)(dest.data() + data_pos));
//        else if (val_type == 2)
//          std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (std::int16_t*)(dest.data() + data_pos));
//        else if (val_type == 3)
//          std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (std::int32_t*)(dest.data() + data_pos));
//        else if (val_type == 4)
//          std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (std::int64_t*)(dest.data() + data_pos));
//        else if (val_type == 5)
//          std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (float*)(dest.data() + data_pos));
////          else if (val_type == 6)
////            std::copy(geno.value_data(), geno.value_data() + geno.non_zero_size(), (double*)(dest.data() + data_pos));
//        else
//          throw std::runtime_error("Unsupported FORMAT type!");
//      }
//    }
  }
}

#endif // LIBSAVVY_BCF_WRITER_HPP