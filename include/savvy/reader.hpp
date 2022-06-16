/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_READER_HPP
#define LIBSAVVY_READER_HPP

#include "site_info.hpp"
#include "savvy.hpp"
#include "file.hpp"
#include "csi.hpp"
#include "s1r.hpp"

#include <shrinkwrap/zstd.hpp>
#include <shrinkwrap/gz.hpp>
#include <shrinkwrap/stdio.hpp>

#include <cstdlib>
#include <string>
#include <memory>
#include <stdexcept>

namespace savvy
{
  //namespace v2
  //{
    class reader : public file
    {
    private:
      std::unique_ptr<std::streambuf> sbuf_;
      std::unique_ptr<std::istream> input_stream_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> ids_;
      typed_value extra_typed_value_;

      std::vector<std::size_t> subset_map_;
      std::size_t subset_size_;

      // Random access
      struct s1r_query_context
      {
        genomic_region reg;
        s1r::reader::query query;
        s1r::reader::query::iterator iter;
        bounding_point bounding_type;
        std::uint32_t current_offset_in_block;
        std::uint32_t total_in_block;
        std::uint64_t total_records_read;
        std::uint64_t max_records_to_read;

        s1r_query_context(s1r::reader& file, genomic_region bounds, bounding_point bound_type = bounding_point::beg) :
          reg(bounds),
          query(file.create_query(bounds)),
          iter(query.begin()),
          bounding_type(bound_type),
          current_offset_in_block(0),
          total_in_block(0),
          total_records_read(0),
          max_records_to_read(std::numeric_limits<std::uint64_t>::max())
        {
        }
      };

      struct csi_query_context
      {
        genomic_region reg;
        std::list<std::pair<std::uint64_t, std::uint64_t>> intervals;
        std::size_t interval_off;
        bounding_point bounding_type;

        csi_query_context(csi_index& file, const std::unordered_map<std::string, std::uint32_t>& contig_map, genomic_region bounds, bounding_point bound_type = bounding_point::beg) :
          reg(bounds),
          interval_off(0),
          bounding_type(bound_type)
        {
          auto tmp = file.query_intervals(reg.chromosome(), contig_map, reg.from(), std::min<std::uint64_t>(reg.to(), std::numeric_limits<std::int64_t>::max())); // TODO: have this return list instead of vector
          intervals.assign(tmp.begin(), tmp.end());
        }
      };

      std::unique_ptr<s1r::reader> s1r_index_;
      std::unique_ptr<s1r_query_context> s1r_query_;
      std::unique_ptr<csi_index> csi_index_;
      std::unique_ptr<csi_query_context> csi_query_;
    public:
      /**
       * Default constuctor.
       */
       reader() : reader("") {}

      /**
       * Constructs reader object and opens SAV, BCF, or VCF file.
       *
       * @param file_path Path to file that will be opened
       */
      reader(const std::string& file_path);

      /**
       * Getter for meta-information lines found in file header.
       *
       * @return Vector of header line key-value pairs
       */
      const std::vector<std::pair<std::string, std::string>>& headers() const { return headers_; }

      /**
       * Getter for sample IDs found in file header
       *
       * @return Vector of sample IDs
       */
      const std::vector<std::string>& samples() const { return ids_; } // TODO: return subset when applicable.

      /**
       * Getter for FORMAT headers found in file.
       *
       * @return List of parsed FORMAT headers
       */
      const std::list<header_value_details>& format_headers() const { return format_headers_; }

      /**
       * Getter for INFO headers found in file.
       *
       * @return List of parsed INFO headers
       */
      const std::list<header_value_details>& info_headers() const { return info_headers_; }

      /**
       * Subsets individual data for future calls to read().
       *
       * @param subset IDs to include if they exist in file
       * @return Ordered vector of file IDs that overlap subset
       */
      std::vector<std::string> subset_samples(const std::unordered_set<std::string>& subset);

      /**
       * Uses S1R or CSI index to query genomic region.
       *
       * @param reg Genomic region to query
       * @param bp Specifies how indels are treated when they cross region bounds
       * @return *this
       */
      reader& reset_bounds(genomic_region reg, bounding_point bp = bounding_point::beg);

      /**
       * Uses S1R index to query records by offset within file.
       *
       * @param reg Bounds of slice query
       * @return *this
       */
      reader& reset_bounds(slice_bounds reg);

      /**
       * Getter for file's phasing status.
       *
       * @return Phasing status
       */
      phasing phasing_status() { return phasing_; }

      /**
       * Sets phasing status. This is useful for ignoring phase when reading BCF/VCF files.
       *
       * @param val Phasing status
       */
      void phasing_status(phasing val) { phasing_ = val; };

      /**
       * Checks for EOF or read error.
       *
       * @return False if either EOF or read error has occurred.
       */
      bool good() const { return this->input_stream_->good(); }

      /**
       * Checks for read error.
       *
       * @return True if read error has occurred
       */
      bool bad() const { return this->input_stream_->bad(); }

      /**
       * Reads next record from file.
       *
       * @param r Destination record object
       * @return *this
       */
      reader& read(variant& r);

      /**
       * Shorthand for read() function.
       *
       * @param r Destination record object
       * @return *this
       */
      reader& operator>>(variant& r) { return read(r); }

      /**
       * Shorthand for good().
       *
       * @return True if good()
       */
      operator bool() const { return good(); };

      /**
       * For SAV files, gets file position for the beginning of current zstd block. For VCF/BCF files, gets "virtual offset".
       *
       * @return File position
       */
      std::streampos tellg() { return this->input_stream_->tellg(); }
    private:
//      void process_header_pair(const std::string& key, const std::string& val);
      bool read_header();
      bool read_header_sav1();

      reader& read_record(variant& r);
      reader& read_vcf_record(variant& r);
      reader& read_sav1_record(variant& r);
      reader& read_indexed_record(variant& r);
      reader& read_csi_indexed_record(variant& r);
    };

    //================================================================//
    // Reader definitions
    inline
    reader::reader(const std::string& file_path)
    {
      FILE* fp = fopen(file_path.c_str(), "rb");
      if (!fp)
      {
        input_stream_ = savvy::detail::make_unique<std::istream>(nullptr);
        return;
      }

      int first_byte = fgetc(fp);
      ungetc(first_byte, fp);

      switch (char(first_byte))
      {
      case '\x1F':
        sbuf_ = ::savvy::detail::make_unique<::shrinkwrap::bgzf::ibuf>(fp);
        break;
      case '\x28':
        sbuf_ = ::savvy::detail::make_unique<::shrinkwrap::zstd::ibuf>(fp);
        break;
      default:
        sbuf_ = ::savvy::detail::make_unique<::shrinkwrap::stdio::filebuf>(fp);
        break;
      }

      input_stream_ = savvy::detail::make_unique<std::istream>(sbuf_.get());

      if (!read_header())
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

      bool csi_exists;
      if (file_format_ == format::sav1 || file_format_ == format::sav2)
        s1r_index_ = ::savvy::detail::make_unique<s1r::reader>(::savvy::detail::file_exists(file_path + ".s1r") ? file_path + ".s1r" : file_path);
      else if ((csi_exists = ::savvy::detail::file_exists(file_path + ".csi")) || ::savvy::detail::file_exists(file_path + ".tbi"))
        csi_index_ = ::savvy::detail::make_unique<csi_index>(file_path + (csi_exists ? ".csi" : ".tbi"));
    }

//    inline
//    sample_subset reader::make_sample_subset(const std::unordered_set<std::string>& subset)
//    {
//      std::vector<std::string> ret;
//      ret.reserve(std::min(subset.size(), ids_.size()));
//
//      std::vector<std::size_t> subset_map(ids_.size(), std::numeric_limits<std::uint64_t>::max());
//
//      std::uint64_t subset_index = 0;
//      for (auto it = ids_.begin(); it != ids_.end(); ++it)
//      {
//        if (subset.find(*it) != subset.end())
//        {
//          subset_map[std::distance(ids_.begin(), it)] = subset_index;
//          ret.push_back(*it);
//          ++subset_index;
//        }
//      }
//
//      return sample_subset{std::move(ret), std::move(subset_map)};
//    }

    inline
    std::vector<std::string> reader::subset_samples(const std::unordered_set<std::string>& subset)
    {
      std::vector<std::string> ret;
      ret.reserve(std::min(subset.size(), ids_.size()));

      subset_map_.clear();
      subset_map_.resize(ids_.size(), std::numeric_limits<std::uint64_t>::max());
      std::uint64_t subset_index = 0;
      for (auto it = ids_.begin(); it != ids_.end(); ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map_[std::distance(ids_.begin(), it)] = subset_index;
          ret.push_back(*it);
          ++subset_index;
        }
      }

      subset_size_ = subset_index;

      return ret;
    }

    inline
    reader& reader::reset_bounds(genomic_region reg, bounding_point bp)
    {
      input_stream_->clear();
      s1r_query_.reset(nullptr);
      csi_query_.reset(nullptr);

      if (s1r_index_ && s1r_index_->good()) //file_format_ == format::sav1 || file_format_ == format::sav2)
      {
        s1r_query_ = ::savvy::detail::make_unique<s1r_query_context>(*s1r_index_, reg, bp);
      }
      else if (csi_index_ && csi_index_->good())
      {
        csi_query_ = ::savvy::detail::make_unique<csi_query_context>(*csi_index_, dict_.str_to_int[dictionary::contig], reg, bp);
        if (!csi_query_->intervals.empty())
        {
          input_stream_->seekg(csi_query_->intervals.front().first);
        }
      }
      else
      {
        input_stream_->setstate(std::ios::failbit); //TODO: error message
      }

      return *this;
    }

    inline
    reader& reader::reset_bounds(slice_bounds reg)
    {
      if (file_format_ == format::sav1 || file_format_ == format::sav2)
      {
        reset_bounds(genomic_region(reg.chromosome()));

        if (s1r_query_) //s1r_index_ && s1r_index_->good())
        {
          auto discard_skip = [this](std::uint32_t num)
          {
            savvy::variant tmp_var;
            while (num > 0 && s1r_query_->current_offset_in_block < s1r_query_->total_in_block && good())
            {
              if (!read_record(tmp_var))
                input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

              ++(s1r_query_->current_offset_in_block);
              --num;
            }
            return num;
          };

          std::size_t num_variants_to_skip = reg.from();
          s1r_query_->max_records_to_read = reg.to() > reg.from() ? reg.to() - reg.from() : 0;
    //        if (num_variants_to_skip < total_in_block_ - current_offset_in_block_)
    //        {
    //          discard_skip(num_variants_to_skip);
    //        }
    //        else
          {
    //          num_variants_to_skip -= (total_in_block_ - current_offset_in_block_);
            while (s1r_query_->iter != s1r_query_->query.end())
            {
              s1r_query_->total_in_block = std::uint32_t(0x000000000000FFFF & s1r_query_->iter->value()) + 1;

              if (num_variants_to_skip < s1r_query_->total_in_block)
              {
                s1r_query_->current_offset_in_block = 0;
                this->input_stream_->seekg(std::streampos((s1r_query_->iter->value() >> 16) & 0x0000FFFFFFFFFFFF));
                ++(s1r_query_->iter);
                discard_skip(num_variants_to_skip);
                return *this;
              }

              num_variants_to_skip -= s1r_query_->total_in_block;
              ++(s1r_query_->iter);
            }

            // Skipped past end of index.
            this->input_stream_->setstate(std::ios::failbit);
          }
        }
      }
      else
      {
        input_stream_->setstate(std::ios::failbit); //TODO: error message
      }

      return *this;
    }

    inline
    reader& reader::read_indexed_record(variant& r)
    {
      while (this->good())
      {
        if (s1r_query_->total_records_read == s1r_query_->max_records_to_read)
        {
          this->input_stream_->setstate(std::ios::eofbit);
          break;
        }

        if (s1r_query_->current_offset_in_block >= s1r_query_->total_in_block)
        {
          if (s1r_query_->iter == s1r_query_->query.end())
          {
            this->input_stream_->setstate(std::ios::eofbit);
            break;
          }
          else
          {
            s1r_query_->total_in_block = std::uint32_t(0x000000000000FFFF & s1r_query_->iter->value()) + 1;
            s1r_query_->current_offset_in_block = 0;
            this->input_stream_->seekg(std::streampos((s1r_query_->iter->value() >> 16) & 0x0000FFFFFFFFFFFF));
            ++(s1r_query_->iter);
          }
        }

        if (!read_record(r))
          input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

        if (!this->good())
        {
          if (s1r_query_->current_offset_in_block < s1r_query_->total_in_block)
          {
            assert(!"Truncated block");
            this->input_stream_->setstate(std::ios::badbit);
          }
        }
        else
        {
          ++(s1r_query_->current_offset_in_block);
          ++(s1r_query_->total_records_read);
          if (region_compare(s1r_query_->bounding_type, r, s1r_query_->reg))
          {
            //this->read_genotypes(annotations, destination);
            break;
          }
          else
          {
            //this->discard_genotypes();
          }
        }
      }
      return *this; //TODO: clear site info before returning if not good
    }

    inline
    reader& reader::read_csi_indexed_record(variant& r)
    {
      while (input_stream_->good())
      {
        if (csi_query_->intervals.empty())
        {
          input_stream_->setstate(std::ios::eofbit);
          break;
        }

        std::uint64_t current_pos = input_stream_->tellg();

        if (current_pos >= csi_query_->intervals.front().second)
        {
          csi_query_->intervals.pop_front();
          if (csi_query_->intervals.empty())
          {
            input_stream_->setstate(std::ios::eofbit);
            return *this;
          }
          else
          {
            assert(csi_query_->intervals.front().first <= csi_query_->intervals.front().second);
            input_stream_->seekg(csi_query_->intervals.front().first);
          }
        }

        //auto pos_before = r.pos();
        if (!read_record(r))
        {
          assert(!"Read failed before end of csi block");
          input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
        }

        //assert(r.pos() >= pos_before);

        if (region_compare(csi_query_->bounding_type, r, csi_query_->reg))
        {
          //this->read_genotypes(annotations, destination);
          break;
        }
        else if (r.pos() > csi_query_->reg.to())
        {
          input_stream_->setstate(std::ios::eofbit);
        }
      }
      return *this; //TODO: clear site info before returning if not good
    }

    inline
    reader& reader::read(variant& r)
    {
      if (good())
      {
        if (s1r_query_)
          return read_indexed_record(r);

        if (csi_query_)
          return read_csi_indexed_record(r);

        if (!read_record(r) && input_stream_->good())
          input_stream_->setstate(std::ios::badbit);
      }

      return *this;
    }

    inline
    reader& reader::read_vcf_record(variant& r)
    {
      if (input_stream_->peek() < 0)
        input_stream_->setstate(input_stream_->rdstate() | std::ios::eofbit);
      else if (!site_info::deserialize_vcf(r, *input_stream_, dict_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
      else if (ids_.size() && !variant::deserialize_vcf2(r, *input_stream_, dict_, ids_.size(), phasing_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
      else
      {
        // TODO: Set not_minimized flag and move minimize routine to writer.
        for (auto it = r.info_.begin(); it != r.info_.end(); ++it)
          it->second.minimize();

        for (auto it = r.format_fields_.begin(); it != r.format_fields_.end(); ++it)
          it->second.minimize();

        if (input_stream_->eof() && (bool)(*input_stream_))
          input_stream_->clear();
      }

      return *this;
    }

    inline
    reader& reader::read_sav1_record(variant& r)
    {
      if (input_stream_->peek() < 0)
        input_stream_->setstate(input_stream_->rdstate() | std::ios::eofbit);
      else if (!site_info::deserialize_sav1(r, *input_stream_, info_headers_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
      else if (!variant::deserialize_sav1(r, *input_stream_, format_headers_, ids_.size()))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
      else
      {
        // TODO: Set not_minimized flag and move minimize routine to writer.
        for (auto it = r.info_.begin(); it != r.info_.end(); ++it)
          it->second.minimize();

        for (auto it = r.format_fields_.begin(); it != r.format_fields_.end(); ++it)
          it->second.minimize();
      }

      return *this;
    }

    inline
    reader& reader::read_record(variant& r)
    {
      if (good())
      {

        if (file_format_ == format::vcf)
          read_vcf_record(r);
        else if (file_format_ == format::sav1)
          read_sav1_record(r);
        else
        {
          std::uint32_t shared_sz, indiv_sz;
          if (!input_stream_->read((char*)&shared_sz, sizeof(shared_sz))) // TODO: set to bad if gcount > 0.
          {
            return *this;
          }

          if (!input_stream_->read((char*)&indiv_sz, sizeof(indiv_sz)))
          {
            std::fprintf(stderr, "Error: Invalid record data\n");
            input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
            return *this;
          }

          if (endianness::is_big())
          {
            shared_sz = endianness::swap(shared_sz);
            indiv_sz = endianness::swap(indiv_sz);
          }

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
          // Read and parse shared and individual data
//          r.shared_data_.resize(shared_sz);
//          r.indiv_buf_.resize(indiv_sz);
//          if (!input_stream_->read(r.shared_data_.data(), r.shared_data_.size()) || !input_stream_->read(r.indiv_buf_.data(), r.indiv_buf_.size()))
//          {
//            std::fprintf(stderr, "Error: Invalid record data\n");
//            input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
//            return *this;
//          }

          std::uint32_t shared_n_samples{};
          if (site_info::deserialize_shared(r, *input_stream_, dict_, shared_n_samples) != shared_sz)
          {
            std::fprintf(stderr, "Error: Invalid shared data\n");
            input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
            return *this;
          }

          bool pbwt_reset = (file_format_ != format::bcf) && (0x800000u & shared_n_samples);
          if (pbwt_reset)
            sort_context_.reset();

          if (variant::deserialize_indiv(r, *input_stream_, dict_, ids_.size(), file_format_ == format::bcf, phasing_) != indiv_sz)
          {
            std::fprintf(stderr, "Error: Invalid individual data\n");
            input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
            return *this;
          }

          if (file_format_ != format::bcf)
          {
            //variant::pbwt_unsort_typed_values(r, extra_typed_value_, sort_context_);
            for (auto it = r.format_fields_.begin(); it != r.format_fields_.end(); ++it)
            {
              if (it->second.pbwt_flag())
              {
                typed_value::internal::delta_decode(it->second, extra_typed_value_, delta_prev_vecs_[it->first]);
                std::swap(it->second, extra_typed_value_);
              }
            }
          }
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        }

        if (good())
        {
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
          // Apply sample subset
          if (subset_size_ != ids_.size()) // TODO: maybe do this after region_compare.
          {
            for (auto it = r.format_fields_.begin(); it != r.format_fields_.end(); ++it)
            {
              it->second.subset(subset_map_, subset_size_, extra_typed_value_);
              it->second.minimize(); // TODO: Set not_minimized flag and move minimize routine to writer.
            }
          }
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        }
      }

      return *this;
    }

    inline
    bool reader::read_header_sav1()
    {
//      std::string version_string(7, '\0');
//      input_stream_->read(&version_string[0], version_string.size());

      input_stream_->read((char*)uuid_.data(), uuid_.size());
      bool parse_ploidy = uuid_.front() != 0;

      std::istreambuf_iterator<char> in_it(*input_stream_);
      std::istreambuf_iterator<char> end;

      std::uint64_t headers_size;
      if (good() && varint_decode(in_it, end, headers_size) != end)
      {
        ++in_it;
        headers_.reserve(1 + headers_size);
        headers_.emplace_back("fileformat","VCFv4.2");

        std::unordered_set<std::string> unique_info_fields;

        while (headers_size && in_it != end)
        {
          std::uint64_t key_size;
          if (varint_decode(in_it, end, key_size) != end)
          {
            ++in_it;
            if (key_size)
            {
              std::string key;
              key.resize(key_size);
              input_stream_->read(&key[0], key_size);

              std::uint64_t val_size;
              if (varint_decode(in_it, end, val_size) != end)
              {
                ++in_it;
                if (key_size)
                {
                  std::string val;
                  val.resize(val_size);
                  input_stream_->read(&val[0], val_size);

                  process_header_pair(key, val);
                  headers_.emplace_back(std::move(key), std::move(val));
                }
              }

            }
          }
          --headers_size;
        }

        if (format_headers_.empty() || (parse_ploidy && std::atoi(format_headers_.back().number.c_str()) == 0))
        {
          return false;
        }

        if (!headers_size)
        {
          std::uint64_t sample_size;
          if (varint_decode(in_it, end, sample_size) != end)
          {
            ++in_it;
            ids_.reserve(sample_size);
            subset_size_ = sample_size;

            std::uint64_t id_sz;
            while (sample_size && varint_decode(in_it, end, id_sz) != end)
            {
              ++in_it;
              ids_.emplace_back();
              if (id_sz)
              {
                ids_.back().resize(id_sz);
                input_stream_->read(&ids_.back()[0], id_sz);
              }
              --sample_size;
            }

            if (!sample_size)
              return true;
          }
        }
      }

      input_stream_->peek();
      return input_stream_->good();
    }

    inline
    bool reader::read_header()
    {
      std::uint32_t header_block_sz = std::uint32_t(-1);

      dict_.str_to_int[dictionary::id]["PASS"] = dict_.entries[dictionary::id].size();
      dict_.entries[dictionary::id].emplace_back(dictionary::entry{"PASS", "", 0});

      std::istream& ifs(*input_stream_);
      int first_byte = ifs.peek();
      if (first_byte == '#')
      {
        file_format_ = format::vcf;
      }
      else
      {
        std::string magic(5, '\0');
        ifs.read(&magic[0], magic.size());
        //transform(magic.begin(), magic.begin() + 3, magic.begin(), ::toupper);
        if (magic[0] == 'B')
        {
          file_format_ = format::bcf;
          ifs.read((char *) &header_block_sz, sizeof(header_block_sz));
        }
        else if (magic.compare(0, 4, "SAV\x02") == 0)
        {
          file_format_ = format::sav2;
          ifs.read((char *) &header_block_sz, sizeof(header_block_sz));
        }
        else if (magic.compare(0, 3, "sav") == 0)
        {
          file_format_ = format::sav1;
          std::array<char, 2> discard;
          ifs.read(discard.data(), 2);
          return read_header_sav1();
        }
        else
        {
          std::fprintf(stderr, "Unsupported file format\n");
          return false;
        }
      }

      if (endianness::is_big())
        header_block_sz = endianness::swap(header_block_sz);

      std::int64_t bytes_read = 0;
      std::string hdr_line;
      while (std::getline(ifs, hdr_line))
      {
        bytes_read += hdr_line.size() + 1;
        if (hdr_line.size() > 1 && hdr_line[1] == 'C')
        {
          std::size_t tab_pos = 0;
          std::size_t last_pos = tab_pos;

          std::int64_t tab_cnt = std::count(hdr_line.begin(), hdr_line.end(), '\t');
          std::int64_t sample_size = std::max<std::int64_t>(tab_cnt - 8, 0);
          ids_.reserve(sample_size);

          if (sample_size)
          {
            while ((tab_pos = hdr_line.find('\t', tab_pos)) != std::string::npos)
            {
              if (tab_cnt < sample_size)
              {
                ids_.emplace_back(hdr_line.substr(last_pos, tab_pos - last_pos));
              }
              last_pos = ++tab_pos;
              --tab_cnt;
            }

            assert(tab_cnt == 0);

            ids_.emplace_back(hdr_line.substr(last_pos, tab_pos - last_pos)); // TODO: allow for no samples.
          }

          assert(ids_.size() == std::size_t(sample_size));
          subset_size_ = sample_size;

          if (header_block_sz - bytes_read < 0)
            break;

          std::array<char, 64> discard;
          while (header_block_sz - bytes_read > 0 && file_format_ != format::vcf)
          {
            auto gcount = ifs.read(discard.data(), std::min(std::size_t(header_block_sz - bytes_read), discard.size())).gcount();
            if (gcount == 0)
              return false;
            bytes_read += gcount;
          }

          return true;
        }

        if (hdr_line.size() < 2)
        {
          break;
        }

        if (hdr_line[1] == '#')
        {
          auto equal_it = std::find(hdr_line.begin(), hdr_line.end(), '=');
          std::string key(hdr_line.begin() + 2, equal_it);
          if (equal_it == hdr_line.end())
          {
            break;
          }
          std::string val(equal_it + 1, hdr_line.end());

          process_header_pair(key, val);
          headers_.emplace_back(std::move(key), std::move(val));
        }
      }

      std::fprintf(stderr, "Error: corrupt header\n");
      return false;
    }

//    inline
//    void reader::process_header_pair(const std::string& key, const std::string& val)
//    {
//      auto hval = parse_header_value(val);
//      if (!hval.id.empty())
//      {
//        int which_dict = -1;
//        if (key == "contig") which_dict = dictionary::contig;
//        else if (key == "INFO" || key == "FILTER" || key == "FORMAT") which_dict = dictionary::id;
//        else if (key == "SAMPLE") which_dict = dictionary::sample;
//
//        if (which_dict >= 0 && dict_.str_to_int[which_dict].find(hval.id) == dict_.str_to_int[which_dict].end())
//        {
//          dictionary::entry e;
//          e.id = hval.id;
//          e.number = hval.number; // TODO: handle special character values.
//          if (hval.type == "Integer")
//            e.type = typed_value::int32;
//          else if (hval.type == "Float")
//            e.type = typed_value::real;
//          else if (hval.type == "String")
//            e.type = typed_value::str;
//
//          if (!hval.idx.empty())
//          {
//            std::size_t idx = std::atoi(hval.idx.c_str());
//            dict_.entries[which_dict].resize(idx, {"DELETED", "", 0});
//          }
//
//          dict_.str_to_int[which_dict][hval.id] = dict_.entries[which_dict].size();
//          dict_.entries[which_dict].emplace_back(std::move(e));
//        }
//      }
//
//      if (key == "INFO")
//      {
//        if (info_headers_map_.find(hval.id) == info_headers_map_.end())
//        {
//          info_headers_.emplace_back(hval);
//          info_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
//        }
//
//        if (hval.id.substr(0, 10) == "_PBWT_SORT")
//        {
//          ::savvy::internal::pbwt_sort_format_context ctx;
//          ctx.format = parse_header_sub_field(val, "Format");
//          ctx.id = hval.id;
//
//          auto insert_it = sort_context_.format_contexts.insert(std::make_pair(std::string(hval.id), std::move(ctx)));
//          sort_context_.field_to_format_contexts.insert(std::make_pair(insert_it.first->second.format, &(insert_it.first->second)));
//        }
//      }
//      else if (key == "FORMAT")
//      {
//        if (format_headers_map_.find(hval.id) == format_headers_map_.end())
//        {
//          format_headers_.emplace_back(hval);
//          format_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
//        }
//      }
//      else if (key == "phasing")
//      {
//        if (val == "none")
//          phasing_ = phasing::none;
//        else if (val == "partial")
//          phasing_ = phasing::partial;
//        else if (val == "phased" || val == "full")
//          phasing_ = phasing::phased;
//      }
//
//      headers_.emplace_back(std::move(key), std::move(val));
//    }
    //================================================================//
  //}

#if 0
  //################################################################//
  class reader_base
  {
  public:
    class sample_iterator
    {
    public:
      typedef sample_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef std::string value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::bidirectional_iterator_tag iterator_category;
    public:
      sample_iterator() : cstring_itr_(nullptr), stdstring_itr_{} {}
      sample_iterator(const char** cstring_itr) : cstring_itr_(cstring_itr), stdstring_itr_{} {}
      sample_iterator(std::vector<std::string>::const_iterator stdstring_itr) : cstring_itr_(nullptr), stdstring_itr_(stdstring_itr) {}

      self_type& operator+=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ += n;
        else
          stdstring_itr_ += n;
        return *this;
      }
      self_type operator+(difference_type n) const { self_type ret(*this); return (ret += n); }
      self_type& operator-=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ -= n;
        else
          stdstring_itr_ -= n;
        return *this;
      }
      self_type operator-(difference_type n) const { self_type ret(*this); return (ret -= n); }
      difference_type operator-(const self_type& b) const
      {
        if (cstring_itr_)
          return cstring_itr_ - b.cstring_itr_;
        return stdstring_itr_ - b.stdstring_itr_;
      }

      self_type& operator--()
      {
        if (cstring_itr_)
          --cstring_itr_;
        else
          --stdstring_itr_;
        return *this;
      }
      self_type operator--(int) { self_type r = *this; --(*this); return r; }
      self_type& operator++()
      {
        if (cstring_itr_)
          ++cstring_itr_;
        else
          ++stdstring_itr_;
        return *this;
      }
      self_type operator++(int) { self_type r = *this; ++(*this); return r; }
      reference operator*()
      {
        if (cstring_itr_)
        {
          tmp_ = *cstring_itr_;
          return tmp_;
        }
        return *stdstring_itr_;
      }
      pointer operator->() { return &(operator*()); }
      bool operator==(const self_type& rhs)
      {
        if (cstring_itr_)
          return cstring_itr_ == rhs.cstring_itr_;
        return stdstring_itr_ == rhs.stdstring_itr_;
      }
      bool operator!=(const self_type& rhs) { return !(*this == rhs); }
    private:
      const char** cstring_itr_;
      std::vector<std::string>::const_iterator stdstring_itr_;
      std::string tmp_;
    };

    virtual ~reader_base() {}

    operator bool() const
    {
      return this->good();
    }

    bool good() const
    {
      if (sav_impl())
        return sav_impl()->good();
      else if (vcf_impl())
        return vcf_impl()->good();
      return false;
    }

    bool fail() const
    {
      if (sav_impl())
        return sav_impl()->fail();
      else if (vcf_impl())
        return vcf_impl()->fail();
      return true;
    }

    bool bad() const
    {
      if (sav_impl())
        return sav_impl()->bad();
      else if (vcf_impl())
        return vcf_impl()->bad();
      return true;
    }

    bool eof() const
    {
      if (sav_impl())
        return sav_impl()->eof();
      else if (vcf_impl())
        return vcf_impl()->eof();
      return true;
    }

//    template <typename T>
//    bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
//    {
//      if (sav_impl())
//        sav_impl()->read_variant(destination, missing_value);
//      else if (vcf_impl())
//        vcf_impl()->read_variant(destination, missing_value);
//
//      return good();
//    }

    const std::vector<std::string>& info_fields() const;
    const std::vector<std::string>& samples() const;
    const std::vector<std::pair<std::string, std::string>>& headers() const;
    std::vector<std::string> subset_samples(const std::set<std::string>& subset);
    void set_policy(enum vcf::empty_vector_policy p);
  private:
    static const std::vector<std::string> empty_string_vector;
    static const std::vector<std::pair<std::string, std::string>> empty_string_pair_vector;
  protected:
    virtual savvy::sav::reader_base* sav_impl() const = 0;
    virtual savvy::vcf::reader_base<1>* vcf_impl() const = 0;
  };
  //################################################################//

  //################################################################//
  class reader : public reader_base
  {
  public:
    reader() {}
    reader(const std::string& file_path, savvy::fmt data_format);
    ~reader() {}

    template <typename T>
    reader& operator>>(variant<T>& destination);

    template <typename T>
    reader& read(site_info& annotations, T& destination);
  private:
    savvy::sav::reader_base* sav_impl() const { return sav_reader_.get(); }
    savvy::vcf::reader_base<1>* vcf_impl() const { return vcf_reader_.get(); }
  private:
    std::unique_ptr<sav::reader> sav_reader_;
    std::unique_ptr<vcf::reader<1>> vcf_reader_;
  };
  //################################################################//

  //################################################################//
  class indexed_reader : public reader_base
  {
  public:
    indexed_reader() {}
    indexed_reader(const std::string& file_path, const genomic_region& reg, savvy::fmt data_format);
    indexed_reader(const std::string& file_path, const genomic_region& reg, bounding_point bounding_type, savvy::fmt data_format);

    void reset_bounds(const genomic_region& reg);
    [[deprecated("Use reset_bounds() instead")]]
    void reset_region(const genomic_region& reg);

    std::vector<std::string> chromosomes() const;

    template <typename T>
    indexed_reader& operator>>(variant<T>& destination);

    template <typename T>
    indexed_reader& read(site_info& annotations, T& destination);

    template <typename Pred, typename T>
    indexed_reader& read_if(Pred fn, site_info& annotations, T& destination);
  private:
    savvy::sav::reader_base* sav_impl() const { return sav_reader_.get(); }
    savvy::vcf::reader_base<1>* vcf_impl() const { return vcf_reader_.get(); }
  private:
    std::unique_ptr<sav::indexed_reader> sav_reader_;
    std::unique_ptr<vcf::indexed_reader<1>> vcf_reader_;
  };
  //################################################################//


  //################################################################//
  template <typename T>
  reader& reader::operator>>(variant<T>& destination)
  {
    return this->read(destination, destination.data());
  }

  template <typename T>
  reader& reader::read(site_info& annotations, T& destination)
  {
    if (sav_impl())
      sav_reader_->read(annotations, destination);
    else if (vcf_impl())
      vcf_reader_->read(annotations, destination);
    return *this;
  }
  //################################################################//

  //################################################################//
  template <typename T>
  indexed_reader& indexed_reader::operator>>(variant<T>& destination)
  {
    return this->read(destination, destination.data());
  }

  template <typename T>
  indexed_reader& indexed_reader::read(site_info& annotations, T& destination)
  {
    if (sav_impl())
      sav_reader_->read(annotations, destination);
    else if (vcf_impl())
      vcf_reader_->read(annotations, destination);
    return *this;
  }

  template <typename Pred, typename T>
  indexed_reader& indexed_reader::read_if(Pred fn, site_info& annoations, T& destination)
  {
    if (sav_reader_)
      sav_reader_->read_if(fn, annoations, destination);
    else if (vcf_reader_)
      vcf_reader_->read_if(fn, annoations, destination);

    return *this;
  }
  //################################################################//
#endif
}

#endif //VC_READER_HPP
