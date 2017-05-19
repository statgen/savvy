#ifndef LIBSAVVY_CMF_READER_HPP
#define LIBSAVVY_CMF_READER_HPP

#include "allele_status.hpp"
#include "varint.hpp"
#include "s1r.hpp"
#include "allele_vector.hpp"
#include "genotype_vector.hpp"
#include "region.hpp"
#include "variant_iterator.hpp"

#include <xzbuf.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <fstream>
#include <tuple>
#include <cmath>
#include <unordered_map>

namespace savvy
{
  namespace sav
  {

    class reader_base
    {
    public:

      reader_base(const std::string& file_path);
      reader_base(reader_base&& source);
      reader_base& operator=(reader_base&& source);
      //reader(const reader&) = delete;
      //reader& operator=(const reader&) = delete;
      virtual ~reader_base() {}

      template <typename T>
      bool read_variant(allele_vector<T>& destination, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN())
      {
        read_variant_details(destination);
        read_genotype(destination, missing_value);

        return good();
      }
      explicit operator bool() const { return input_stream_.good(); }
      bool good() const { return input_stream_.good(); }
      bool fail() const { return input_stream_.fail(); }
      bool bad() const { return input_stream_.bad(); }
      std::uint64_t sample_count() const { return this->sample_ids_.size(); }
      std::uint64_t haplotype_count() const { return this->sample_count() * this->ploidy(); }
      std::vector<std::string>::const_iterator samples_begin() const { return sample_ids_.begin(); }
      std::vector<std::string>::const_iterator samples_end() const { return sample_ids_.end(); }
//      std::vector<std::string>::const_iterator prop_fields_begin() const { return metadata_fields_.begin(); }
//      std::vector<std::string>::const_iterator prop_fields_end() const { return metadata_fields_.end(); }
      const std::string& chromosome() const { return chromosome_; }
      std::vector<std::string> chromosomes() const { return {chromosome_}; }
      std::vector<std::string> prop_fields() const { return std::vector<std::string>(metadata_fields_); }
      std::uint8_t ploidy() const { return ploidy_level_; }
      const std::string& file_path() const { return file_path_; }
      std::streampos tellg() { return this->input_stream_.tellg(); }
    protected:
      virtual bool update_file_position() { return true; }
      template <typename T>
      void read_variant_details(allele_vector<T>& destination)
      {
        if (good())
        {
          if (!this->update_file_position())
          {
            this->input_stream_.setstate(std::ios::failbit);
          }
          else
          {
            std::istreambuf_iterator<char> in_it(input_stream_);
            std::istreambuf_iterator<char> end_it;

            std::uint64_t locus;
            if (varint_decode(in_it, end_it, locus) == end_it)
            {
              this->input_stream_.setstate(std::ios::badbit);
            }
            else
            {
              ++in_it;
              std::uint64_t sz;
              if (varint_decode(in_it, end_it, sz) == end_it)
              {
                this->input_stream_.setstate(std::ios::badbit);
              }
              else
              {
                ++in_it;
                std::string ref;
                ref.resize(sz);
                if (sz)
                  input_stream_.read(&ref[0], sz);

                if (varint_decode(in_it, end_it, sz) == end_it)
                {
                  this->input_stream_.setstate(std::ios::badbit);
                }
                else
                {
                  ++in_it;
                  std::string alt;
                  alt.resize(sz);
                  if (sz)
                    input_stream_.read(&alt[0], sz);

                  std::unordered_map<std::string, std::string> props;
                  props.reserve(this->metadata_fields_.size());
                  std::string prop_val;
                  for (const std::string& key : metadata_fields_)
                  {
                    if (varint_decode(in_it, end_it, sz) == end_it)
                    {
                      this->input_stream_.setstate(std::ios::badbit);
                      break;
                    }
                    else
                    {
                      ++in_it;
                      if (sz)
                      {
                        prop_val.resize(sz);
                        input_stream_.read(&prop_val[0], sz);
                        props[key] = prop_val;
                      }
                    }
                  }

                  destination = allele_vector<T>(std::string(chromosome()), locus, std::move(ref), std::move(alt), std::move(props), std::move(destination));
                  destination.resize(0);
                }
              }
            }
          }
        }
      }

      template <typename T>
      void read_genotype(allele_vector<T>& destination, const typename T::value_type missing_value)
      {
        if (good())
        {
          const typename T::value_type alt_value = typename T::value_type(1);
          std::istreambuf_iterator<char> in_it(input_stream_);
          std::istreambuf_iterator<char> end_it;

          destination.resize(sample_count() * ploidy());

          std::uint64_t sz;
          varint_decode(in_it, end_it, sz);
          std::uint64_t total_offset = 0;
          for (std::size_t i = 0; i < sz && in_it != end_it; ++i, ++total_offset)
          {
            std::uint8_t allele;
            std::uint64_t offset;
            one_bit_prefixed_varint::decode(++in_it, end_it, allele, offset);
            total_offset += offset;
            destination[total_offset] = (allele ? missing_value : alt_value);
          }

          input_stream_.get();
        }
      }
    protected:
      std::vector<std::string> sample_ids_;
      std::string chromosome_;
      //ixzbuf sbuf_;
      ixzstream input_stream_;
      std::string file_path_;
      std::uint8_t ploidy_level_;
      std::vector<std::string> metadata_fields_;
    };

    class reader : public reader_base
    {
    public:
      using reader_base::reader_base;

      template <typename T>
      reader& operator>>(allele_vector<T>& destination)
      {
        read_variant(destination);
        return *this;
      }
      template <typename T>
      reader& read(allele_vector<T>& destination, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN())
      {
        read_variant(destination, missing_value);
        return *this;
      }
    };

    class indexed_reader : public reader_base
    {
    public:
      indexed_reader(const std::string& file_path, std::uint64_t from, std::uint64_t to, const std::string& index_file_path = "") :
        reader_base(file_path),
        index_(index_file_path.size() ? index_file_path : file_path + ".s1r"),
        query_(index_.create_query(from, to)),
        i_(query_.begin())
      {
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      indexed_reader(const std::string& file_path, const region& reg, const std::string& index_file_path = "") :
        indexed_reader(file_path, reg.from(), reg.to(), index_file_path)
      {
        if (this->chromosome() != reg.chromosome())
          this->input_stream_.setstate(std::ios::badbit);
      }

      template <typename T>
      indexed_reader& operator>>(allele_vector<T>& destination)
      {
        read_variant(destination);
        return *this;
      }

      template <typename T>
      indexed_reader& read(allele_vector<T>& destination, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN())
      {
        read_variant(destination, missing_value);
        return *this;
      }

      template <typename T, typename Pred>
      indexed_reader& read_if(allele_vector<T>& destination, Pred fn, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN())
      {
        bool predicate_failed = true;
        while (good() && predicate_failed)
        {
          read_variant_details(destination);
          if (good())
          {
            predicate_failed = !fn(destination);
            if (!predicate_failed)
            {
              read_genotype(destination, missing_value);
            }
          }
        }

        return *this;
      }

      void reset_region(std::uint64_t from, std::uint64_t to)
      {
        input_stream_.clear();
        query_ = index_.create_query(from, to);
        i_ = query_.begin();
        if (!index_.good())
          this->input_stream_.setstate(std::ios::badbit);
      }

      void reset_region(const region& reg)
      {
        input_stream_.clear();
        reset_region(reg.from(), reg.to());
        if (this->chromosome() != reg.chromosome())
          this->input_stream_.setstate(std::ios::badbit);
      }

    private:
      bool update_file_position()
      {
        if (i_ == query_.end())
          return false;
        input_stream_.seekg(std::streampos(i_->value().first));
        ++i_;
        return true;
      }
    private:
      s1r::reader index_;
      s1r::reader::query query_;
      s1r::reader::query::iterator i_;
    };

    class writer
    {
    public:
      template <typename RandAccessStringIterator, typename RandAccessStringIterator2>
      writer(const std::string& file_path, const std::string& chromosome, std::uint8_t ploidy, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end, RandAccessStringIterator2 property_fields_beg, RandAccessStringIterator2 property_fields_end) :
        output_stream_(file_path),
        file_path_(file_path),
        sample_size_(samples_end - samples_beg),
        ploidy_level_(ploidy),
        record_count_(0),
        block_size_(64)
      {
        std::string version_string("sav\x00\x01\x00\x00", 7);
        output_stream_.write(version_string.data(), version_string.size());

        std::ostreambuf_iterator<char> out_it(output_stream_);

        varint_encode(chromosome.size(), out_it);
        std::copy(chromosome.begin(), chromosome.end(), out_it);
        varint_encode(ploidy_level_, out_it);

        varint_encode(sample_size_, out_it);
        for (auto it = samples_beg; it != samples_end; ++it)
        {
          std::size_t str_sz = get_string_size(*it);
          varint_encode(str_sz, out_it);
          if (str_sz)
            output_stream_.write(&(*it)[0], str_sz);
        }

        property_fields_.assign(property_fields_beg, property_fields_end);
        varint_encode(property_fields_.size(), out_it);
        for (auto it = property_fields_.begin(); it != property_fields_.end(); ++it)
        {
          std::size_t str_sz = get_string_size(*it);
          varint_encode(str_sz, out_it);
          if (str_sz)
            output_stream_.write(&(*it)[0], str_sz);
        }
      }


      template <typename RandAccessStringIterator>
      writer(const std::string& file_path, const std::string& chromosome, std::uint8_t ploidy, RandAccessStringIterator samples_beg, RandAccessStringIterator samples_end) :
        writer(file_path, chromosome, ploidy, std::forward<RandAccessStringIterator>(samples_beg), std::forward<RandAccessStringIterator>(samples_end), empty_string_array.end(), empty_string_array.end())
      {

      }

      template <typename T>
      writer& operator<<(const allele_vector<T>& m)
      {
        if (output_stream_.good())
        {
          if (m.size() != sample_size_ * ploidy_level_)
          {
            output_stream_.setstate(std::ios::failbit);
          }
          else
          {
            if ((record_count_ % block_size_) == 0)
              output_stream_.flush();
            write<T>(m);
            ++record_count_;
          }
        }
        return *this;
      }

//#define NO_LEB128 1
#ifdef NO_LEB128
      template <typename T>
      void write(const allele_vector<T>& m)
      {
        const typename T::value_type ref_value = typename T::value_type();
        //std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
        std::uint64_t sz = m.locus();

        output_stream_.write((char*)&sz, 8); //varint_encode(m.locus(), os_it);

        sz = (m.ref().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          output_stream_.write(m.ref().data(), m.ref().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        sz = (m.alt().size() << 48);
        output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
        if (m.alt().size())
          output_stream_.write(m.alt().data(), m.alt().size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        for (const std::string& key : property_fields_)
        {
          std::string value(m.prop(key));
          sz = (m.ref().size() << 48);
          output_stream_.write((char*)&sz, 2); //varint_encode(m.ref().size(), os_it);
          if (value.size())
            output_stream_.write(value.data(), value.size()); //std::copy(m.ref().begin(), m.ref().end(), os_it);
          //os.write(&source.ref_[0], source.ref_.size());
        }

        struct sparse_geno
        {
          std::uint32_t v: 1, offset: 31;
        };

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        output_stream_.write((char*)&non_zero_count, 8);//varint_encode(non_zero_count, os_it);

        std::vector<sparse_geno> tmp(non_zero_count);

        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        std::size_t non_ref_counter = 0;
        for (auto it = beg; it != m.end(); ++it)
        {
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            tmp[non_ref_counter].v = (std::isnan(*it)  ? std::uint8_t(0x80) : std::uint8_t(0x00));
            tmp[non_ref_counter].offset = offset;
            ++non_ref_counter;
          }
        }
        output_stream_.write((char*)tmp.data(), tmp.size() * 4);
      }
#else
      template <typename T>
      void write(const allele_vector<T>& m)
      {
        const typename T::value_type ref_value = typename T::value_type();
        std::ostreambuf_iterator<char> os_it(output_stream_.rdbuf());
        varint_encode(m.locus(), os_it);

        varint_encode(m.ref().size(), os_it);
        if (m.ref().size())
          std::copy(m.ref().begin(), m.ref().end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

        varint_encode(m.alt().size(), os_it);
        if (m.alt().size())
          std::copy(m.alt().begin(), m.alt().end(), os_it);
        //os.write(&source.alt_[0], source.alt_.size());

        for (const std::string& key : property_fields_)
        {
          std::string value(m.prop(key));
          varint_encode(value.size(), os_it);
          if (value.size())
            std::copy(value.begin(), value.end(), os_it);
        }

        std::uint64_t non_zero_count =  m.size() - static_cast<std::size_t>(std::count(m.begin(), m.end(), ref_value));
        varint_encode(non_zero_count, os_it);
        std::uint64_t last_pos = 0;
        auto beg = m.begin();
        for (auto it = beg; it != m.end(); ++it)
        {
          if (*it != ref_value)
          {
            std::uint64_t dist = static_cast<std::uint64_t>(std::distance(beg, it));
            std::uint64_t offset = dist - last_pos;
            last_pos = dist + 1;
            std::uint8_t allele = (std::isnan(*it)  ? std::uint8_t(0x80) : std::uint8_t(0x00));
            one_bit_prefixed_varint::encode(allele, offset, os_it);
          }
        }
      }
#endif

      static bool create_index(const std::string& input_file_path, std::string output_file_path = "");
    private:
      template <typename T>
      static std::size_t get_string_size(T str);
    private:
      static const std::array<std::string, 0> empty_string_array;
    private:
      std::ofstream output_stream_;
      std::vector<std::string> property_fields_;
      std::string file_path_;
      std::uint64_t sample_size_;
      std::uint8_t ploidy_level_;
      std::uint32_t metadata_fields_cnt_;
      std::size_t record_count_;
      std::size_t block_size_;
    };

    template <typename T>
    std::size_t writer::get_string_size(T str)
    {
      return str.size();
    }

    template <>
    inline std::size_t writer::get_string_size<const char*>(const char* str)
    {
      return std::strlen(str);
    }

    template <typename VecType>
    using variant_iterator =  basic_variant_iterator<reader_base, VecType>;

    template <typename ValType>
    using dense_variant_iterator =  basic_variant_iterator<reader_base, std::vector<ValType>>;
    template <typename ValType>
    using sparse_variant_iterator =  basic_variant_iterator<reader_base, compressed_vector<ValType>>;
  }
}

#endif //LIBSAVVY_CMF_READER_HPP
