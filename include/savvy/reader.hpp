/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_READER_HPP
#define LIBSAVVY_READER_HPP

#include "site_info.hpp"
#include "sav_reader.hpp"
#include "vcf_reader.hpp"
#include "savvy.hpp"
#include "file.hpp"

#include <shrinkwrap/zstd.hpp>
#include <shrinkwrap/gz.hpp>

#include <cstdlib>
#include <string>
#include <memory>
#include <stdexcept>

namespace savvy
{
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

  namespace v2
  {
    class reader : public file
    {
    private:
      std::string file_path_;
      std::unique_ptr<std::streambuf> sbuf_;
      std::unique_ptr<std::istream> input_stream_;
      std::array<std::uint8_t, 16> uuid_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<header_value_details> info_headers_;
      std::unordered_map<std::string, std::reference_wrapper<header_value_details>> info_headers_map_;
      std::vector<header_value_details> format_headers_;
      std::unordered_map<std::string, std::reference_wrapper<header_value_details>> format_headers_map_;

      std::vector<std::string> ids_;
      dictionary dict_;
      format file_format_;
      phasing phasing_;

      ::savvy::internal::pbwt_sort_context sort_context_;

      // Random access
      struct index_data
      {
        s1r::reader file;
        genomic_region reg;
        s1r::reader::query query;
        s1r::reader::query::iterator iter;
        bounding_point bounding_type;
        std::uint32_t current_offset_in_block;
        std::uint32_t total_in_block;
        std::uint64_t total_records_read;
        std::uint64_t max_records_to_read;

        index_data(const std::string& file_path, genomic_region bounds, bounding_point bound_type = bounding_point::beg) :
          file(file_path),
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

      std::unique_ptr<index_data> index_;
    public:
      reader(const std::string& file_path);

      const std::vector<std::pair<std::string, std::string>>& headers() const { return headers_; }
      const std::vector<std::string>& samples() const { return ids_; }

      const std::vector<header_value_details>& format_headers() const { return format_headers_; }
      const std::vector<header_value_details>& info_headers() const { return info_headers_; }

      /**
       *
       * @param subset IDs to include if they exist in file.
       * @return subset context to be used with typed_value::get().
       */
      sample_subset make_sample_subset(const std::unordered_set<std::string>& subset);

      reader& reset_bounds(genomic_region reg);
      phasing phasing_status() { return phasing_; }

      bool good() const { return this->input_stream_->good(); }
      bool bad() const { return this->input_stream_->bad(); }

      reader& read(variant& r);
      reader& operator>>(variant& r) { return read(r); }

      operator bool() const { return good(); };
    private:
      void process_header_pair(const std::string& key, const std::string& val);
      bool read_header();
      bool read_header_sav1();

      reader& read_record(variant& r);
      reader& read_vcf_record(variant& r);
      reader& read_sav1_record(variant& r);
      reader& read_indexed_record(variant& r);
    };

    //================================================================//
    // Reader definitions
    inline
    reader::reader(const std::string& file_path) :
      file_path_(file_path)
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
        throw std::runtime_error("uncompressed files not yet supported.");
      }

      input_stream_ = savvy::detail::make_unique<std::istream>(sbuf_.get());

      if (!read_header())
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
    }

    inline
    sample_subset reader::make_sample_subset(const std::unordered_set<std::string>& subset)
    {
      std::vector<std::string> ret;
      ret.reserve(std::min(subset.size(), ids_.size()));

      std::vector<std::size_t> subset_map(ids_.size(), std::numeric_limits<std::uint64_t>::max());

      std::uint64_t subset_index = 0;
      for (auto it = ids_.begin(); it != ids_.end(); ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map[std::distance(ids_.begin(), it)] = subset_index;
          ret.push_back(*it);
          ++subset_index;
        }
      }

      return sample_subset{std::move(ret), std::move(subset_map)};
    }

    inline
    reader& reader::reset_bounds(genomic_region reg)
    {

      input_stream_->clear();

      index_ = ::savvy::detail::make_unique<index_data>(::savvy::detail::file_exists(file_path_ + ".s1r") ? file_path_ + ".s1r" : file_path_, reg);
      if (!index_->file.good())
      {
        input_stream_->setstate(input_stream_->rdstate() | std::ios::failbit);
      }

      return *this;
    }

    inline
    reader& reader::read_indexed_record(variant& r)
    {
      while (this->good())
      {
        if (index_->total_records_read == index_->max_records_to_read)
        {
          this->input_stream_->setstate(std::ios::eofbit);
          break;
        }

        if (index_->current_offset_in_block >= index_->total_in_block)
        {
          if (index_->iter == index_->query.end())
            this->input_stream_->setstate(std::ios::eofbit);
          else
          {
            index_->total_in_block = std::uint32_t(0x000000000000FFFF & index_->iter->value()) + 1;
            index_->current_offset_in_block = 0;
            this->input_stream_->seekg(std::streampos((index_->iter->value() >> 16) & 0x0000FFFFFFFFFFFF));
            ++(index_->iter);
          }
        }

        if (!read_record(r))
          input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

        if (!this->good())
        {
          if (index_->current_offset_in_block < index_->total_in_block)
          {
            assert(!"Truncated block");
            this->input_stream_->setstate(std::ios::badbit);
          }
        }
        else
        {
          ++(index_->current_offset_in_block);
          ++(index_->total_records_read);
          if (region_compare(index_->bounding_type, r, index_->reg))
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
    reader& reader::read(variant& r)
    {
      if (good())
      {
        if (index_)
          return read_indexed_record(r);

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
      else if (!variant::deserialize_vcf(r, *input_stream_, dict_, ids_.size(), phasing_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

      return *this;
    }

    inline
    reader& reader::read_sav1_record(variant& r)
    {
      if (input_stream_->peek() < 0)
        input_stream_->setstate(input_stream_->rdstate() | std::ios::eofbit);
      else if (!site_info::deserialize_sav1(r, *input_stream_, info_headers_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
      else if (!variant::deserialize_sav1(r, *input_stream_, format_headers_, ids_.size(), phasing_))
        input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);

      return *this;
    }

    inline
    reader& reader::read_record(variant& r)
    {
      if (good())
      {
        if (file_format_ == format::vcf)
          return read_vcf_record(r);
        else if (file_format_ == format::sav1)
          return read_sav1_record(r);

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

        // TODO: endianess

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Read and parse shared and individual data
        r.shared_data_.resize(shared_sz);
        r.indiv_buf_.resize(indiv_sz);
        if (!input_stream_->read(r.shared_data_.data(), r.shared_data_.size()) || !input_stream_->read(r.indiv_buf_.data(), r.indiv_buf_.size()))
        {
          std::fprintf(stderr, "Error: Invalid record data\n");
          input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
          return *this;
        }

        if (!variant::deserialize(r, dict_, this->ids_.size(), file_format_ == format::bcf, phasing_))
        {
          std::fprintf(stderr, "Error: Invalid record data\n");
          input_stream_->setstate(input_stream_->rdstate() | std::ios::badbit);
          return *this;
        }
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Handle semantic INFO fields
        int pbwt_reset{}; r.get_info("_PBWT_RESET", pbwt_reset);
        if (pbwt_reset)
          sort_context_.reset();

        std::vector<::savvy::internal::pbwt_sort_format_context*> pbwt_format_pointers;
        pbwt_format_pointers.reserve(r.format_fields_.size());
        for (auto it = r.info().begin(); it != r.info().end(); )
        {
          if (it->first.substr(0, 10) == "_PBWT_SORT")
          {
            auto f = sort_context_.format_contexts.find(it->first);
            if (f != sort_context_.format_contexts.end())
            {
              pbwt_format_pointers.emplace_back(&(f->second));
              it = r.remove_info(it);
              continue;
            }
          }
          ++it;
        }
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        // Unsort FMT fields
        for (auto fmt_it = r.format_fields_.begin(); fmt_it != r.format_fields_.end(); ++fmt_it)
        {
          for (auto pt = pbwt_format_pointers.begin(); pt != pbwt_format_pointers.end(); )
          {
            if ((*pt)->format == fmt_it->first)
            {
              typed_value::internal::pbwt_unsort(fmt_it->second, (*pt)->sort_map, sort_context_.prev_sort_mapping, sort_context_.counts);
              pt = pbwt_format_pointers.erase(pt);
              break;
            }
            else
            {
              ++pt;
            }
          }
        }
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

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
        headers_.reserve(headers_size);

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
          std::size_t sample_size = tab_cnt - 8;
          ids_.reserve(sample_size);

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

          assert(ids_.size() == sample_size);

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
        }
      }

      std::fprintf(stderr, "Error: corrupt header\n");
      return false;
    }

    inline
    void reader::process_header_pair(const std::string& key, const std::string& val)
    {
      auto hval = parse_header_value(val);
      if (!hval.id.empty())
      {
        int which_dict = -1;
        if (key == "contig") which_dict = dictionary::contig;
        else if (key == "INFO" || key == "FILTER" || key == "FORMAT") which_dict = dictionary::id;
        else if (key == "SAMPLE") which_dict = dictionary::sample;

        if (which_dict >= 0 && dict_.str_to_int[which_dict].find(hval.id) == dict_.str_to_int[which_dict].end())
        {
          dictionary::entry e;
          e.id = hval.id;
          e.number = hval.number; // TODO: handle special character values.
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

      if (key == "INFO")
      {
        if (info_headers_map_.find(hval.id) == info_headers_map_.end())
        {
          info_headers_.emplace_back(hval);
          info_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
        }

        if (hval.id.substr(0, 10) == "_PBWT_SORT")
        {
          ::savvy::internal::pbwt_sort_format_context ctx;
          ctx.format = parse_header_sub_field(val, "Format");
          ctx.id = hval.id;

          auto insert_it = sort_context_.format_contexts.insert(std::make_pair(std::string(hval.id), std::move(ctx)));
          sort_context_.field_to_format_contexts.insert(std::make_pair(insert_it.first->second.format, &(insert_it.first->second)));
        }
      }
      else if (key == "FORMAT")
      {
        if (format_headers_map_.find(hval.id) == format_headers_map_.end())
        {
          format_headers_.emplace_back(hval);
          format_headers_map_.insert(std::make_pair(hval.id, std::ref(info_headers_.back())));
        }
      }
      else if (key == "phasing")
      {
        if (val == "none")
          phasing_ = phasing::none;
        else if (val == "partial")
          phasing_ = phasing::partial;
        else if (val == "phased" || val == "full")
          phasing_ = phasing::phased;
      }

      headers_.emplace_back(std::move(key), std::move(val));
    }
    //================================================================//
  }
}

#endif //VC_READER_HPP
