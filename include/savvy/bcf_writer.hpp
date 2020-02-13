#ifndef LIBSAVVY_BCF_WRITER_HPP
#define LIBSAVVY_BCF_WRITER_HPP

#include "savvy/utility.hpp"

#include <shrinkwrap/zstd.hpp>

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <chrono>

namespace savvy
{
  namespace bcf
  {
    template<typename T>
    typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, std::uint8_t>::type
    int_type(T val)
    { // TODO: Handle missing values
      if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
        return 0x01;
      else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
        return 0x02;
      else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
        return 0x03;
      else
        return 0x04;
    }

    template <typename Iter, typename IntT>
    Iter deserialize_int(Iter it, Iter end, IntT& dest)
    {
      if (it != end)
      {
        std::uint8_t type_byte = *(it++);
        std::size_t int_width = 1u << bcf_type_shift[type_byte & 0x0Fu];
        if (end - it >= int_width)
        {
          switch (int_width)
          {
          case 1u:
            {
              dest = IntT(*(it++));
              return it;
            }
          case 2u:
            {
              std::int16_t tmp;
              char *tmp_p = (char *)&tmp;
              *(tmp_p) = *(it++);
              *(tmp_p + 1) = *(it++);
              dest = tmp;
              return it;
            }
          case 4u:
            {
              std::int32_t tmp;
              char *tmp_p = (char *)&tmp;
              *(tmp_p) = *(it++);
              *(tmp_p + 1) = *(it++);
              *(tmp_p + 2) = *(it++);
              *(tmp_p + 3) = *(it++);
              dest = tmp;
              return it;
            }
          case 8u:
            {
              std::int64_t tmp;
              char *tmp_p = (char *)&tmp;
              *(tmp_p) = *(it++);
              *(tmp_p + 1) = *(it++);
              *(tmp_p + 2) = *(it++);
              *(tmp_p + 3) = *(it++);
              *(tmp_p + 4) = *(it++);
              *(tmp_p + 5) = *(it++);
              *(tmp_p + 6) = *(it++);
              *(tmp_p + 7) = *(it++);
              dest = tmp;
              return it;
            }
          }
        }
      }
      throw std::runtime_error("Not a BCF integer");
    }

//    template <typename Iter>
//    Iter deserialize_string(Iter it, Iter end, std::string& dest)
//    {
//      if (it == end) return end;
//
//      std::uint8_t type_byte = *(it++);
//      if (it == end || (type_byte & 0x0Fu) != 0x07u)
//        throw std::runtime_error("Not a BCF string");
//
//      std::int32_t sz = (type_byte >> 4u);
//      if (sz == 15)
//        it = deserialize_int(it, end, sz);
//
//      if (end - it < sz)
//        throw std::runtime_error("Invalid byte sequence");
//
//      dest.resize(sz);
//      std::copy_n(it, sz, dest.begin());
//      return it + sz;
//    }

    template <typename Iter, typename VecT>
    typename std::enable_if<std::is_same<typename std::iterator_traits<Iter>::value_type, char>::value, Iter>::type
    deserialize_vec(Iter it, Iter end, VecT& dest)
    {
      if (it == end)
        throw std::runtime_error("Invalid byte sequence");

      std::uint8_t type_byte = *(it++);

      std::int32_t sz = (type_byte >> 4u);
      if (sz == 15)
        it = deserialize_int(it, end, sz);

      std::size_t type_width = 1u << bcf_type_shift[0x0Fu & type_byte];

      if (end - it < sz * type_width)
        throw std::runtime_error("Invalid byte sequence");

      dest.resize(sz);
      char* char_p = &(*it);
      switch (0x0Fu & type_byte)
      {
      case 0x01u:
      {
        auto p = (std::int8_t *)char_p;
        std::copy_n(p, sz, dest.begin());
      }
      case 0x02u:
      {
        auto p = (std::int16_t *)char_p;
        std::copy_n(p, sz, dest.begin()); // TODO: use transform to switch endian when needed
      }
      case 0x03u:
      {
        auto p = (std::int32_t *)char_p;
        std::copy_n(p, sz, dest.begin());
      }
      case 0x04u:
      {
        auto p = (std::int64_t *)char_p;
        std::copy_n(p, sz, dest.begin());
      }
      case 0x05u:
      {
        auto p = (float *)char_p;
        std::copy_n(p, sz, dest.begin());
      }
      case 0x07u:
      {
        std::copy_n(char_p, sz, dest.begin());
      }
      }

      return it + (sz * type_width);
    }

    template <typename OutT, typename T>
    typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
    serialize_typed_int_exact(OutT out_it, const T& val)
    {
      T v = static_cast<T>(val);
      std::uint8_t type;
      if (std::is_same<T, std::int8_t>::value)
        type = 0x01;
      else if (std::is_same<T, std::int16_t>::value)
        type = 0x02;
      else if (std::is_same<T, std::int32_t>::value)
        type = 0x03;
      else if (std::is_same<T, std::int64_t>::value)
        type = 0x04;
      else
        return false;

      *(out_it++) = (1u << 4u) | type;

      char* p_end = ((char*)&v) + sizeof(T);
      for (char* p = (char*)&v; p != p_end; ++p)
      {
        *(out_it++) = *p;
      }
      return true;
    }

    template <typename OutT, typename T>
    typename std::enable_if<std::is_signed<T>::value, bool>::type
    serialize_typed_scalar(OutT out_it, const T& val)
    {
      std::uint8_t type_byte = 1;
      if (std::is_integral<T>::value)
      {
        if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
          return serialize_typed_int_exact(out_it, std::int8_t(val));
        else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
          return serialize_typed_int_exact(out_it, std::int16_t(val));
        else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
          return serialize_typed_int_exact(out_it, std::int32_t(val));
        else
          return serialize_typed_int_exact(out_it, std::int64_t(val));
      }

      if (std::is_same<T, float>::value)
        *(out_it++) = (1u << 4u) | 0x05u;
      else if (std::is_same<T, double>::value)
        *(out_it++) = (1u << 4u) | 0x06u;
      else
        return false;

      char* p_end = ((char*)&val) + sizeof(T);
      for (char* p = (char*)&val; p != p_end; ++p)
      {
        *(out_it++) = *p;
      }

      return true;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, bool>::type
    write_typed_scalar(std::ostream& os, const T& val)
    {
      std::uint8_t type_byte = 1;
      if (std::is_integral<T>::value)
      {
        if (std::is_same<T, std::int8_t>::value)
        {
          type_byte = (type_byte << 4) | 0x01;
        }
        else if (std::is_same<T, std::int16_t>::value)
        {
          type_byte = (type_byte << 4) | 0x02;
        }
        else if (std::is_same<T, std::int32_t>::value)
        {
          type_byte = (type_byte << 4) | 0x03;
        }
        else
        {
          return false;
        }
      }
      else if (std::is_same<T, float>::value)
      {
        type_byte = (type_byte << 4) | 0x05;
      }
      else
      {
        return false;
      }

      os.write((char*)&type_byte, 1);
      os.write((char*)&val, sizeof(val));
      return true;
    }

    template <typename OutT>
    bool serialize_type_and_size(OutT out_it, std::uint8_t type, std::size_t size)
    {
      if (size < 15)
      {
        *out_it = size << 4u | type;
        ++out_it;
        return true;
      }

      *out_it = 0xF0 | type;
      ++out_it;

      return serialize_typed_scalar(out_it, (std::int64_t)size);
    }

    template <typename Iter, typename T>
    typename std::enable_if<std::is_signed<T>::value, void>::type
    serialize_typed_vec(Iter out_it, const std::vector<T>& vec) // TODO: use smallest int type.
    {
      static_assert(!std::is_same<T, std::int64_t>::value && !std::is_same<T, double>::value, "64-bit integers not allowed in BCF spec.");

      std::uint8_t type_byte = vec.size() < 15 ? vec.size() : 15;
      if (std::is_same<T, std::int8_t>::value)
      {
        type_byte = (type_byte << 4) | 0x01;
      }
      else if (std::is_same<T, std::int16_t>::value)
      {
        type_byte = (type_byte << 4) | 0x02;
      }
      else if (std::is_same<T, std::int32_t>::value)
      {
        type_byte = (type_byte << 4) | 0x03;
      }
      else if (std::is_same<T, float>::value)
      {
        type_byte = (type_byte << 4) | 0x05;
      }

      *out_it = type_byte;

      if (vec.size() >= 15)
      {
        if (vec.size() <= 0x7F)
          serialize_typed_int_exact(out_it, (std::int8_t)vec.size());
        else if (vec.size() <= 0x7FFF)
          serialize_typed_int_exact(out_it, (std::int16_t)vec.size());
        else if (vec.size() <= 0x7FFFFFFF)
          serialize_typed_int_exact(out_it, (std::int32_t)vec.size());
        else
          throw std::runtime_error("string too big");
      }

      std::copy_n((char*)vec.data(), sizeof(T) * vec.size(), out_it);
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, void>::type
    write_typed_vec(std::ostream& os, const std::vector<T>& vec)
    {
      static_assert(!std::is_same<T, std::int64_t>::value && !std::is_same<T, double>::value, "64-bit integers not allowed in BCF spec.");

      std::uint8_t type_byte = vec.size() < 15 ? vec.size() : 15;
      if (std::is_same<T, std::int8_t>::value)
      {
        type_byte = (type_byte << 4) | 0x01;
      }
      else if (std::is_same<T, std::int16_t>::value)
      {
        type_byte = (type_byte << 4) | 0x02;
      }
      else if (std::is_same<T, std::int32_t>::value)
      {
        type_byte = (type_byte << 4) | 0x03;
      }
      else if (std::is_same<T, float>::value)
      {
        type_byte = (type_byte << 4) | 0x05;
      }

      os.write((char*)&type_byte, 1);

      if (vec.size() >= 15)
      {
        if (vec.size() <= 0x7F)
          write_typed_scalar(os, (std::int8_t)vec.size());
        else if (vec.size() <= 0x7FFF)
          write_typed_scalar(os, (std::int16_t)vec.size());
        else if (vec.size() <= 0x7FFFFFFF)
          write_typed_scalar(os, (std::int32_t)vec.size());
        else
          throw std::runtime_error("vector too big");
      }

      os.write((char*)vec.data(), std::int32_t(sizeof(T) * vec.size()));
    }

    template <typename OutT>
    void serialize_typed_str(OutT out_it, const std::string& str)
    {
      std::uint8_t type_byte = str.size() < 15 ? str.size() : 15;
      type_byte = (type_byte << 4) | 0x07;

      *out_it = type_byte;

      if (str.size() >= 15)
      {
        if (str.size() <= 0x7F)
          serialize_typed_int_exact(out_it, (std::int8_t)str.size());
        else if (str.size() <= 0x7FFF)
          serialize_typed_int_exact(out_it, (std::int16_t)str.size());
        else if (str.size() <= 0x7FFFFFFF)
          serialize_typed_int_exact(out_it, (std::int32_t)str.size());
        else
          throw std::runtime_error("string too big");
      }

      std::copy_n(str.begin(), str.size(), out_it);
    }

    inline void write_typed_str(std::ostream& os, const std::string& str)
    {
      std::uint8_t type_byte = str.size() < 15 ? str.size() : 15;
      type_byte = (type_byte << 4) | 0x07;


      os.write((char*)&type_byte, 1);
      if (str.size() >= 15)
      {
        if (str.size() <= 0x7F)
          serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int8_t)str.size());
        else if (str.size() <= 0x7FFF)
          serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int16_t)str.size());
        else if (str.size() <= 0x7FFFFFFF)
          serialize_typed_int_exact(std::ostreambuf_iterator<char>(os), (std::int32_t)str.size());
        else
          throw std::runtime_error("string too big");
      }

      os.write(str.data(), str.size());
    }

    template <typename T>
    typename std::enable_if<std::is_signed<typename T::value_type>::value, std::uint32_t>::type
    get_typed_value_size(const T& vec)
    {
      static_assert(!std::is_same<typename T::value_type, std::int64_t>::value && !std::is_same<typename T::value_type, double>::value, "64-bit integers not allowed in BCF spec.");
      std::uint32_t ret;
      if (vec.size() < 15)
        ret = 1;
      else if (vec.size() <= 0x7F)
        ret = 2 + 1;
      else if (vec.size() <= 0x7FFF)
        ret = 2 + 2;
      else if (vec.size() <= 0x7FFFFFFF)
        ret = 2 + 4;
      else
        return -1; // vec too big

        ret += vec.size() * sizeof(typename T::value_type);
      return ret;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint32_t>::type
    get_typed_value_size(T vec)
    {
      static_assert(!std::is_same<T, std::int64_t>::value, "64-bit integers not allowed in BCF spec.");
      return 1 + sizeof(T);
    }


//    struct record
//    {
//      std::uint32_t l_shared = 0;
//      std::uint32_t l_indiv = 0;
//      std::int32_t chrom_idx = 1;
//      std::int32_t pos = 123;
//      std::int32_t ref_len = 1;
//      float qual = 1;
//      std::uint32_t n_allele_info = (2 << 16) | 0;
//      std::uint32_t n_fmt_sample = (1 << 24) | 6;
//      std::string variant_id = "rs123";
//      std::vector<std::string> alleles = {"C","T"};
//      std::vector<std::string> filters;
//      std::vector<std::pair<std::string, std::string>> info;
//    };
//
//    class writer
//    {
//    private:
//      std::string file_path_;
//      std::ofstream file_;
//      std::vector<std::pair<std::string, std::string>> headers_;
//      std::vector<std::string> samples_;
//    public:
//      writer(const std::string& file_path,
//        std::vector<std::pair<std::string, std::string>> headers,
//        std::vector<std::string> samples) :
//        file_path_(file_path),
//        file_(file_path, std::ios::binary),
//        headers_(std::move(headers)),
//        samples_(std::move(samples))
//      {
//
//      }
//
//      void write_header()
//      {
//        std::string magic = "BCF\x02\x02";
//
//        std::uint32_t header_block_sz = 0;
//        for (auto it = headers_.begin(); it != headers_.end(); ++it)
//        {
//          header_block_sz += it->first.size();
//          header_block_sz += it->second.size();
//          header_block_sz += 4;
//        }
//
//        std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
//        header_block_sz += column_names.size();
//
//        std::vector<std::string> samples = {"1","2","3","4","5","6"};
//        for (auto it = samples.begin(); it != samples.end(); ++it)
//        {
//          header_block_sz += 1 + it->size();
//        }
//
//        header_block_sz += 2; //new line and null
//
//        file_.write(magic.data(), magic.size());
//        file_.write((char*)(&header_block_sz), sizeof(header_block_sz));
//
//        for (auto it = headers_.begin(); it != headers_.end(); ++it)
//        {
//          file_.write("##", 2);
//          file_.write(it->first.data(), it->first.size());
//          file_.write("=", 1);
//          file_.write(it->second.data(), it->second.size());
//          file_.write("\n", 1);
//        }
//
//        file_.write(column_names.data(), column_names.size());
//        for (auto it = samples.begin(); it != samples.end(); ++it)
//        {
//          file_.write("\t", 1);
//          file_.write(it->data(), it->size());
//        }
//
//        file_.write("\n\0", 2);
//
//      }
//
//      writer& write(record r, std::vector<std::int8_t>& genos)
//      {
//        r.l_shared += sizeof(r.chrom_idx);
//        r.l_shared += sizeof(r.pos);
//        r.l_shared += sizeof(r.ref_len);
//        r.l_shared += sizeof(r.qual);
//        r.l_shared += sizeof(r.n_allele_info);
//        r.l_shared += sizeof(r.n_fmt_sample);
//        r.l_shared += get_typed_value_size(r.variant_id);
//
//        for (auto& a : r.alleles)
//          r.l_shared += get_typed_value_size(a);
//
//        std::vector<std::int8_t> filter_vec = {0};
//        r.l_shared += get_typed_value_size(filter_vec);
//        for (auto& i: r.info)
//          r.l_shared += (get_typed_value_size(0) + get_typed_value_size(i.second));
//
//        r.l_indiv += get_typed_value_size(std::int8_t(0));
//        r.l_indiv += 1;
//        r.l_indiv += genos.size() * sizeof(std::vector<std::int8_t>::value_type);
//
//        file_.write((char*)(&r.l_shared), sizeof(r.l_shared));
//        file_.write((char*)(&r.l_indiv), sizeof(r.l_indiv));
//        file_.write((char*)(&r.chrom_idx), sizeof(r.chrom_idx));
//        file_.write((char*)(&r.pos), sizeof(r.pos));
//        file_.write((char*)(&r.ref_len), sizeof(r.ref_len));
//        file_.write((char*)(&r.qual), sizeof(r.qual));
//        file_.write((char*)(&r.n_allele_info), sizeof(r.n_allele_info));
//        file_.write((char*)(&r.n_fmt_sample), sizeof(r.n_fmt_sample));
//
//
//        write_typed_str(file_, r.variant_id);
//        write_typed_str(file_, r.alleles[0]);
//        write_typed_str(file_, r.alleles[1]);
//        write_typed_vec(file_, filter_vec);
//        write_typed_scalar(file_, std::int8_t(3));
//
//        std::uint16_t ploidy = 2;
//        std::uint8_t typed_byte = (ploidy << 4u) | 1u;
//        file_.write((char*)&typed_byte, 1);
//        file_.write((char*)genos.data(), genos.size());
//
//        return *this;
//      }
//    };
//
//    inline void write_bcf()
//    {
//      std::ofstream ofs("test.bcf", std::ios::binary);
//
//      // BEGIN HEADER
//      std::string magic = "BCF\x02\x02";
//
//      std::vector<std::pair<std::string, std::string>> headers = {
//        {"fileformat","VCFv4.3"},
//        {"FILTER","<ID=PASS,Description=\"Pass\">"},
//        {"FILTER","<ID=q10,Description=\"Quality below 10\">"},
//        {"FILTER","<ID=s50,Description=\"Less than 50% of samples have data\">"},
//        {"FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">"},
//        {"fileDate","20090805"},
//        {"source","myImputationProgramV3.1"},
//        {"reference","file:///seq/references/1000GenomesPilot-NCBI36.fasta"},
//        {"contig","<ID=18,length=80373285,assembly=B36,species=\"Homo sapiens\",taxonomy=x>"},
//        {"contig","<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"},
//        {"phasing","partial"},
//        {"INFO","<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"},
//        {"INFO","<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"},
//        {"INFO","<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"},
//        {"INFO","<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"},
//        {"INFO","<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"},
//        {"INFO","<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"}};
//
//
//      std::uint32_t header_block_sz = 0;
//      for (auto it = headers.begin(); it != headers.end(); ++it)
//      {
//        header_block_sz += it->first.size();
//        header_block_sz += it->second.size();
//        header_block_sz += 4;
//      }
//
//      std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
//      header_block_sz += column_names.size();
//
//      std::vector<std::string> samples = {"1","2","3","4","5","6"};
//      for (auto it = samples.begin(); it != samples.end(); ++it)
//      {
//        header_block_sz += 1 + it->size();
//      }
//
//      header_block_sz += 2; //new line and null
//
//      ofs.write(magic.data(), magic.size());
//      ofs.write((char*)(&header_block_sz), sizeof(header_block_sz));
//
//      for (auto it = headers.begin(); it != headers.end(); ++it)
//      {
//        ofs.write("##", 2);
//        ofs.write(it->first.data(), it->first.size());
//        ofs.write("=", 1);
//        ofs.write(it->second.data(), it->second.size());
//        ofs.write("\n", 1);
//      }
//
//      ofs.write(column_names.data(), column_names.size());
//      for (auto it = samples.begin(); it != samples.end(); ++it)
//      {
//        ofs.write("\t", 1);
//        ofs.write(it->data(), it->size());
//      }
//
//      ofs.write("\n\0", 2);
//      // END HEADER
//
//      // BEGIN RECORDS
//      record r;
//      r.l_shared += sizeof(r.chrom_idx);
//      r.l_shared += sizeof(r.pos);
//      r.l_shared += sizeof(r.ref_len);
//      r.l_shared += sizeof(r.qual);
//      r.l_shared += sizeof(r.n_allele_info);
//      r.l_shared += sizeof(r.n_fmt_sample);
//      r.l_shared += get_typed_value_size(r.variant_id);
//      for (auto& a : r.alleles)
//        r.l_shared += get_typed_value_size(a);
//      std::vector<std::int8_t> filter_vec = {0};
//      r.l_shared += get_typed_value_size(filter_vec);
//      for (auto& i: r.info)
//        r.l_shared += (get_typed_value_size(0) + get_typed_value_size(i.second));
//
//      r.l_indiv += get_typed_value_size(std::int8_t(0));
//      r.l_indiv += 1;
//      r.l_indiv += 6 * 2 * sizeof(char);
//
//      ofs.write((char*)(&r.l_shared), sizeof(r.l_shared));
//      ofs.write((char*)(&r.l_indiv), sizeof(r.l_indiv));
//      ofs.write((char*)(&r.chrom_idx), sizeof(r.chrom_idx));
//      ofs.write((char*)(&r.pos), sizeof(r.pos));
//      ofs.write((char*)(&r.ref_len), sizeof(r.ref_len));
//      ofs.write((char*)(&r.qual), sizeof(r.qual));
//      ofs.write((char*)(&r.n_allele_info), sizeof(r.n_allele_info));
//      ofs.write((char*)(&r.n_fmt_sample), sizeof(r.n_fmt_sample));
//
//
//      write_typed_str(ofs, r.variant_id);
//      write_typed_str(ofs, r.alleles[0]);
//      write_typed_str(ofs, r.alleles[1]);
//      write_typed_vec(ofs, filter_vec);
//      write_typed_scalar(ofs, std::int8_t(3));
//
//      std::uint8_t typed_byte = (2 << 4) | 1;
//      ofs.write((char*)&typed_byte, 1);
//      ofs.write("\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02", 12);
//
//
//
//
////      std::int64_t scalar64_test;
////      get_typed_value_size(scalar64_test);
//      // END RECORDS
//
//    }
  }

  namespace sav2
  {
    class typed_value
    {
    public:
      static const std::uint8_t int8   = 1;
      static const std::uint8_t int16  = 2;
      static const std::uint8_t int32  = 3;
      static const std::uint8_t int64  = 4;
      static const std::uint8_t real   = 5;
      static const std::uint8_t real64 = 6;
      static const std::uint8_t str    = 7;
      static const std::uint8_t sparse = 8;

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
      is_missing(const T& v);

      template <typename T>
      static typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value, bool>::type
      is_missing(const T& v);

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type
      type_code();

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type
      type_code(const T& val);

      template <typename T>
      static typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type
      type_code_ignore_missing(const T& val);

      enum class get_status : std::uint8_t
      {
        ok = 0,
        does_not_fit,
        not_a_scalar,
        not_a_vector
      };
    public:
      typed_value() { }

      template<typename T>
      typed_value(const T& v)
      {
        init(v);
      }

      typed_value(std::uint8_t type, std::size_t sz, char* data_ptr);

      typed_value(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char* data_ptr);

      typed_value(typed_value&& src)
      {
        operator=(std::move(src));
      }

      typed_value& operator=(typed_value&& src);

      template <typename T>
      bool operator>>(savvy::compressed_vector<T>& dest) const
      {
        if (this->off_ptr_)
        {
        }
        else if (this->val_ptr_)
        {

        }
        return false;
      }

      class internal
      {
      public:
        template <typename Iter>
        static void serialize(const typed_value& v, Iter out_it);
      };
    private:
      template<typename T>
      typename std::enable_if<std::is_signed<T>::value, void>::type
      init(const T& v);

      template <typename T>
      struct is_dense_vector
      {
        static const bool value =
          std::is_signed<typename T::value_type>::value
            && !std::is_same<std::string, T>::value
            &&
            (
              std::is_same<std::random_access_iterator_tag, typename std::iterator_traits<typename T::iterator>::iterator_category>::value
#if 0
            || std::is_same<T, boost::numeric::ublas::vector<typename T::value_type>>::value
#endif
          );
      };

      template<typename T>
      typename std::enable_if<is_dense_vector<T>::value, void>::type
      init(const T& vec);

      template<typename T>
      typename std::enable_if<std::is_same<T, ::savvy::compressed_vector<typename T::value_type>>::value && std::is_signed<typename T::value_type>::value, void>::type
      init(const T& vec);

      template<typename T>
      typename std::enable_if<std::is_same<T, std::string>::value, void>::type
      init(const T& vec);
    private:
      std::uint8_t val_type_ = 0;
      std::uint8_t off_type_ = 0;
      std::size_t size_ = 0;
      std::size_t sparse_size_ = 0;
      char* val_ptr_ = nullptr;
      char* off_ptr_ = nullptr;
      std::vector<char> local_data_;
    };

    class dictionary
    {
    public:
      static const std::uint8_t id = 0;
      static const std::uint8_t contig = 1;
      static const std::uint8_t sample = 2;
      std::array<std::unordered_map<std::string, std::uint32_t>, 3> str_to_int;
      std::array<std::vector<std::string>, 3> int_to_str;
    };

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

    class site_info
    {
    private:
      std::string chrom_;
      std::string id_;
      std::uint32_t pos_ = (std::uint32_t)-1;
      float qual_ = 0.f;
      std::string ref_;
      std::vector<std::string> alts_;
      std::vector<std::string> filters_;
      std::vector<std::pair<std::string, typed_value>> info_;
      std::vector<char> shared_data_;
    protected:
      std::uint32_t n_fmt_ = 0;
    public:
      site_info() {}
      site_info(std::string chrom, std::uint32_t pos, std::string ref, std::vector<std::string> alts,
        std::string id = "",
        float qual = 0.f,
        std::vector<std::string> filters = {},
        std::vector<std::pair<std::string, typed_value>> info = {});

      const std::string& chrom() const { return chrom_; }
      const std::string& id() const { return id_; }
      std::uint32_t pos() const { return pos_ + 1u; }
      float qual() const { return qual_; }
      const std::string& ref() const { return ref_; }
      const std::vector<std::string>& alts() const { return alts_; }
      const std::vector<std::string>& filters() const { return filters_; }
      const std::vector<std::pair<std::string, typed_value>>& info() const { return info_; }

      class internal
      {
      public:
        static bool read(site_info& s, std::istream& ifs, const dictionary& dict, std::uint32_t shared_sz);

        template <typename Itr>
        static bool serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint16_t n_fmt);
      };
    };

    class variant : public site_info
    {
    public:
      class internal
      {
      public:
        static bool read(variant& v, std::istream& ifs, const dictionary& dict, std::size_t sample_size, bool is_bcf);

        static bool serialize(const variant& v, std::vector<char>& buf, const dictionary& dict, std::size_t sample_size, bool is_bcf, std::uint32_t& shared_sz, std::uint32_t& indiv_sz);
      };
    private:
      std::vector<std::pair<std::string, typed_value>> format_fields_;
      std::vector<char> indiv_buf_;
    public:
      using site_info::site_info;
      const decltype(format_fields_)& format_fields() const { return format_fields_; }

      template <typename T>
      void set_format(const std::string& key, const std::vector<T>& geno, std::set<std::string> sparse_keys = {"GT","EC","DS","HDS"});

      template <typename T>
      void set_format(const std::string& key, const compressed_vector<T>& geno);
    };


    class reader
    {
    private:
      std::ifstream ifs_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> ids_;
      dictionary dict_;
    public:
      reader(const std::string& file_path) :
        ifs_(file_path, std::ios::binary)
      {
        if (!read_header(ifs_, headers_, ids_, dict_))
          ifs_.setstate(ifs_.rdstate() | std::ios::badbit);
      }

      reader& read_record(variant& r)
      {
        bool is_bcf = false; // TODO ...
        if (!variant::internal::read(r, ifs_, dict_, ids_.size(), is_bcf))
          ifs_.setstate(ifs_.rdstate() | std::ios::badbit);

        return *this;
      }

      operator bool() const { return (bool)ifs_; };
    private:
      static bool read_header(std::istream& ifs, std::vector<std::pair<std::string, std::string>>& headers, std::vector<std::string>& ids, dictionary& dict)
      {
        std::string magic(5, '\0');
        ifs.read(&magic[0], magic.size());

        std::uint32_t header_block_sz = 0;
        ifs.read((char*)&header_block_sz, sizeof(header_block_sz));

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
            ids.reserve(sample_size);

            while ((tab_pos = hdr_line.find('\t', tab_pos)) != std::string::npos)
            {
              if (tab_cnt < sample_size)
              {
                ids.emplace_back(hdr_line.substr(last_pos, tab_pos - last_pos));
              }
              last_pos = ++tab_pos;
              --tab_cnt;
            }

            assert(tab_cnt == 0);

            ids.emplace_back(hdr_line.substr(last_pos, tab_pos - last_pos));

            assert(ids.size() == sample_size);

            if (header_block_sz - bytes_read < 0)
              break;

            std::array<char, 64> discard;
            while (header_block_sz - bytes_read > 0)
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

            std::string hid = parse_header_sub_field(val, "ID");
            if (!hid.empty())
            {
              int which_dict = key == "contig" ? dictionary::contig : dictionary::id;

              dict.int_to_str[which_dict].emplace_back(hid);
              dict.str_to_int[which_dict][hid] = dict.int_to_str[which_dict].size() - 1;
            }

            headers.emplace_back(std::move(key), std::move(val));
          }
        }

        std::fprintf(stderr, "Error: corrupt header");
        return false;
      }
    };


    class writer
    {
    private:
      std::mt19937_64 rng_;
      std::unique_ptr<std::streambuf> output_buf_;
      std::string file_path_;
      std::ostream ofs_;
      std::array<std::uint8_t, 16> uuid_;
      dictionary dict_;
      std::size_t n_samples_ = 0;
      std::vector<char> serialized_buf_;

      // Data members to support indexing
      std::unique_ptr<s1r::writer> index_file_;
      std::string current_chromosome_;
      std::size_t block_size_ = 4096;
      std::size_t record_count_ = 0;
      std::size_t record_count_in_block_ = 0;
      std::uint32_t current_block_min_ = std::numeric_limits<std::uint32_t>::max();
      std::uint32_t current_block_max_ = 0;
    private:
      static std::filebuf* create_std_filebuf(const std::string& file_path, std::ios::openmode mode)
      {
        std::filebuf* ret = new std::filebuf();
        ret->open(file_path.c_str(), mode);
        return ret;
      }

      static std::unique_ptr<std::streambuf> create_out_streambuf(const std::string& file_path, std::int8_t compression_level)
      {
        if (compression_level > 0)
          return std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path, compression_level));
        else
          return std::unique_ptr<std::streambuf>(create_std_filebuf(file_path, std::ios::binary | std::ios::out));
      }
    public:
      writer(const std::string& file_path) :
        rng_(std::chrono::high_resolution_clock::now().time_since_epoch().count() ^ std::clock() ^ (std::uint64_t)this),
        output_buf_(create_out_streambuf(file_path, 3)), //opts.compression == compression_type::zstd ? std::unique_ptr<std::streambuf>(new shrinkwrap::zstd::obuf(file_path)) : std::unique_ptr<std::streambuf>(new std::filebuf(file_path, std::ios::binary))),
        ofs_(output_buf_.get()),
        //samples_(samples_beg, samples_end),
        file_path_(file_path),
        uuid_(::savvy::detail::gen_uuid(rng_))
      {
        // TODO: Use mkstemp when shrinkwrap supports FILE*
        if (true && file_path_ != "/dev/stdout") // TODO: Check if indexing and zstd is enabled.
        {
          index_file_ = ::savvy::detail::make_unique<s1r::writer>(file_path + ".s1r" , uuid_);
        }
      }
      void write_header(const std::vector<std::pair<std::string, std::string>>& headers, const std::vector<std::string>& ids)
      {
        std::string magic = {'S','A','V','\x02','\x00'};

        std::uint32_t header_block_sz = 0;
        for (auto it = headers.begin(); it != headers.end(); ++it)
        {
          header_block_sz += it->first.size();
          header_block_sz += it->second.size();
          header_block_sz += 4;
        }

        std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        header_block_sz += column_names.size();

        for (auto it = ids.begin(); it != ids.end(); ++it)
        {
          header_block_sz += 1 + it->size();
        }

        header_block_sz += 2; //new line and null

        ofs_.write(magic.data(), magic.size());
        ofs_.write((char*)(&header_block_sz), sizeof(header_block_sz));

        for (auto it = headers.begin(); it != headers.end(); ++it)
        {
          std::string hid = parse_header_sub_field(it->second, "ID");
          if (!hid.empty())
          {
            int which_dict = it->first == "contig" ? dictionary::contig : dictionary::id;

            dict_.int_to_str[which_dict].emplace_back(hid);
            dict_.str_to_int[which_dict][hid] = dict_.int_to_str[which_dict].size() - 1;
          }

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

        ofs_.write("\n\0", 2);
      }

      writer& write_record(const variant& r)
      {
        bool is_bcf = true; // TODO: ...

        if (block_size_ != 0 && ((record_count_ % block_size_) == 0 || r.chrom() != current_chromosome_))
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
        }

        std::uint32_t shared_sz, indiv_sz;
        if (!variant::internal::serialize(r, serialized_buf_, dict_, n_samples_, is_bcf, shared_sz, indiv_sz))
        {
          ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
        }
        else
        {
          ofs_.write((char*)&shared_sz, sizeof(shared_sz));
          ofs_.write((char*)&indiv_sz, sizeof(indiv_sz));
          ofs_.write(serialized_buf_.data(), serialized_buf_.size());

          current_block_min_ = std::min(current_block_min_, std::uint32_t(r.pos()));
          std::size_t max_alt_size = 0;
          for (auto it = r.alts().begin(); it != r.alts().end(); ++it)
            max_alt_size = std::max(max_alt_size, it->size());
          current_block_max_ = std::max(current_block_max_, std::uint32_t(r.pos() + std::max(r.ref().size(), max_alt_size)) - 1);
          ++record_count_in_block_;
          ++record_count_;
        }

        return *this;

        //variant::internal::write(r, ofs_, );
//
//        std::uint32_t l_shared = 6 * 4; // chrom through n.fmt.sample
//
//        l_shared += bcf::get_typed_value_size(r.id());
//
//        l_shared += bcf::get_typed_value_size(r.ref());
//        for (auto& a : r.alts())
//          l_shared += bcf::get_typed_value_size(a);
//
//        std::vector<std::int8_t> filter_vec(r.filters().size()); // TODO: Allow more than 255 filters.
//        for (std::size_t i = 0; i < filter_vec.size(); ++i)
//        {
//          auto it = dict_.str_to_int[dictionary::id].find(r.filters()[i]);
//          if (it == dict_.str_to_int[dictionary::id].end())
//          {
//            std::cerr << "Error: filter not valid" << std::endl;
//            ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
//          }
//          assert(it->second <= 0xFF);
//          filter_vec[i] = it->second;
//        }
//
//        l_shared += bcf::get_typed_value_size(filter_vec);
//        for (auto& i: r.info())
//        {
//          auto it = dict_.str_to_int[dictionary::contig].find(i.first);
//          if (it == dict_.str_to_int[dictionary::contig].end())
//          {
//            std::cerr << "Error: INFO key not valid" << std::endl;
//            ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
//            return;
//          }
//
//          l_shared += bcf::get_typed_value_size((std::int32_t)it->second); // TODO: Allow for variable sized INFO keys.
//          l_shared += i.second.data.size(); //(bcf::get_typed_value_size(0) + bcf::get_typed_value_size(i.second));
//        }
//
//
//        std::uint32_t l_indiv = 0;
//        for (auto& f : r.format_fields())
//        {
//          auto it = dict_.str_to_int[dictionary::id].find(f.key());
//          if (it != dict_.str_to_int[dictionary::id].end())
//          {
//            l_indiv += bcf::get_typed_value_size((std::int32_t)it->second);
//            l_indiv += f.serialized_data().size();
//          }
//        }
//
//        ofs_.write((char*)(&l_shared), sizeof(l_shared));
//        ofs_.write((char*)(&l_indiv), sizeof(l_indiv));
//
//        {
//          auto it = dict_.str_to_int[dictionary::contig].find(r.chrom());
//          if (it == dict_.str_to_int[dictionary::contig].end())
//          {
//            std::cerr << "Error: contig not valid" << std::endl;
//            ofs_.setstate(ofs_.rdstate() | std::ios::badbit);
//            return;
//          }
//          ofs_.write((char*)(&it->second), sizeof(it->second));
//        }
//
//
//        std::uint32_t tmp_i32 = r.pos();
//        ofs_.write((char*)(&tmp_i32), sizeof(tmp_i32));
//        tmp_i32 = r.ref().size();
//        ofs_.write((char*)(&tmp_i32), sizeof(tmp_i32));
//        float qual = r.qual();
//        ofs_.write((char*)(&qual), sizeof(qual));
//
//        tmp_i32 = ((1 + r.alts().size()) << 16) | r.info().size();
//        ofs_.write((char*)(&tmp_i32), sizeof(tmp_i32));
//        tmp_i32 = (r.format_fields().size() << 24) | 0; // zero for n_sample in SAV
//        ofs_.write((char*)(&tmp_i32), sizeof(tmp_i32));
//
//
//        bcf::write_typed_str(ofs_, r.id());
//        bcf::write_typed_str(ofs_, r.ref());
//        for (auto& a: r.alts())
//          bcf::write_typed_str(ofs_, a);
//        bcf::write_typed_vec(ofs_, filter_vec);
//
//        std::int16_t i = 0;
//        float f = static_cast<float>(i);
//        i = static_cast<std::uint16_t>(f);
//
//        for (auto& i: r.info())
//        {
//          bcf::write_typed_scalar(ofs_, (std::int32_t)dict_.str_to_int[dictionary::id].find(i.first)->second); // TODO: Allow for variable sized INFO keys.
//          bcf::serialize_type_and_size(std::ostreambuf_iterator<char>(ofs_), i.second.type, i.second.length);
//          ofs_.write((char*)(i.second.data.data()), i.second.data.size());
//        }
//
//        for (auto& f : r.format_fields())
//        {
//          bcf::write_typed_scalar(ofs_, (std::int32_t)(dict_.str_to_int[dictionary::id].find(f.key())->second)); // TODO: Allow for variable sized FMT keys.
//          ofs_.write((char*)(f.serialized_data().data()), f.serialized_data().size());
//        }
      }
    };

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value && std::is_integral<T>::value, bool>::type
    typed_value::is_missing(const T& v)
    {
      return v == std::numeric_limits<T>::min();
    }

    template <typename T>
    typename std::enable_if<std::is_same<T, float>::value || std::is_same<T, double>::value, bool>::type
    typed_value::is_missing(const T& v)
    {
      return std::isnan(v);
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type typed_value::type_code()
    {
      if (std::is_same<T, std::int8_t>::value)
        return typed_value::int8;
      if (std::is_same<T, std::int16_t>::value)
        return typed_value::int16;
      if (std::is_same<T, std::int32_t>::value)
        return typed_value::int32;
      if (std::is_same<T, std::int64_t>::value)
        return typed_value::int64;
      if (std::is_same<T, float>::value)
        return typed_value::real;
      if (std::is_same<T, double>::value)
        return typed_value::real64;
      return 0;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type typed_value::type_code_ignore_missing(const T& val)
    {
      std::uint8_t type = type_code<T>();
      if (type >= typed_value::int16 && type <= typed_value::int64)
      {
        if (val <= std::numeric_limits<std::int8_t>::max() && val >= std::numeric_limits<std::int8_t>::min())
          type = typed_value::int8;
        else if (val <= std::numeric_limits<std::int16_t>::max() && val >= std::numeric_limits<std::int16_t>::min())
          return typed_value::int16;
        else if (val <= std::numeric_limits<std::int32_t>::max() && val >= std::numeric_limits<std::int32_t>::min())
          return typed_value::int32;
        else
          return typed_value::int64;
      }
      return type;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint8_t>::type typed_value::type_code(const T& val)
    {
      std::uint8_t type = type_code<T>();
      if (type >= typed_value::int16 && type <= typed_value::int64)
      {
        if (val <= std::numeric_limits<std::int8_t>::max() && val > std::numeric_limits<std::int8_t>::min())
          type = typed_value::int8;
        else if (val <= std::numeric_limits<std::int16_t>::max() && val > std::numeric_limits<std::int16_t>::min())
          return typed_value::int16;
        else if (val <= std::numeric_limits<std::int32_t>::max() && val > std::numeric_limits<std::int32_t>::min())
          return typed_value::int32;
        else
          return typed_value::int64;
      }
      return type;
    }

    inline
    typed_value::typed_value(std::uint8_t type, std::size_t sz, char* data_ptr) :
      val_type_(type),
      size_(sz),
      val_ptr_(data_ptr)
    {
    }

    inline
    typed_value::typed_value(std::uint8_t val_type, std::size_t sz, std::uint8_t off_type, std::size_t sp_sz, char* data_ptr) :
      val_type_(val_type),
      off_type_(off_type),
      size_(sz),
      sparse_size_(sp_sz),
      val_ptr_(data_ptr + sp_sz * (1u << bcf_type_shift[off_type])),
      off_ptr_(data_ptr)
    {
    }

    inline
    typed_value& typed_value::operator=(typed_value&& src)
    {
      if (&src != this)
      {
        val_type_ = src.val_type_;
        off_type_ = src.off_type_;
        size_ = src.size_;
        sparse_size_ = src.sparse_size_;
        val_ptr_ = src.val_ptr_;
        off_ptr_ = src.off_ptr_;
        local_data_ = std::move(src.local_data_);

        src.val_type_ = 0;
        src.off_type_ = 0;
        src.size_ = 0;
        src.sparse_size_ = 0;
        src.val_ptr_ = nullptr;
        src.off_ptr_ = nullptr;
      }
      return *this;
    }

    template <typename Iter>
    void typed_value::internal::serialize(const typed_value& v, Iter out_it)
    {

      std::uint8_t type_byte =  v.off_type_ ? typed_value::sparse : v.val_type_;
      type_byte = std::uint8_t(std::min(std::size_t(15), v.size_) << 4u) | type_byte;
      *(out_it++) = type_byte;
      if (v.size_ >= 15u)
        bcf::serialize_typed_scalar(out_it, static_cast<std::int64_t>(v.size_));

      if (v.off_type_)
      {
        type_byte = std::uint8_t(v.off_type_ << 4u) | v.val_type_;
        *(out_it++) = type_byte;
        bcf::serialize_typed_scalar(out_it, static_cast<std::int64_t>(v.sparse_size_));
        std::size_t pair_width = (1u << bcf_type_shift[v.off_type_]) + (1u << bcf_type_shift[v.val_type_]);
        std::copy_n(v.val_ptr_, v.sparse_size_ * pair_width, out_it);
      }
      else
      {
        std::copy_n(v.val_ptr_, v.size_ * (1u << bcf_type_shift[v.val_type_]), out_it);
      }
    }

    template<typename T>
    typename std::enable_if<std::is_signed<T>::value, void>::type
    typed_value::init(const T& v)
    {
      val_type_ = type_code<T>();
      size_ = 1;
      std::size_t width = 1u << bcf_type_shift[val_type_];
      // TODO: handle endianess
      local_data_.resize(width);
      std::memcpy(local_data_.data(), &v, width);
      val_ptr_ = local_data_.data();
    }

    template<typename T>
    typename std::enable_if<typed_value::is_dense_vector<T>::value, void>::type
    typed_value::init(const T& vec)
    {
      typedef typename T::value_type vtype;
      if (std::is_integral<vtype>::value)
      {
        vtype min_val = 0;
        vtype max_val = 0;
        for (auto it = vec.begin(); it != vec.end(); ++it)
        {
          if (!is_missing(*it))
          {
            if (*it > max_val)
              max_val = *it;
            if (*it < min_val)
              min_val = *it;
          }
        }

        val_type_ = type_code(std::min(vtype(-max_val), min_val));
      }
      else
      {
        val_type_ = type_code<vtype>();
      }

      size_ = vec.size();
      std::size_t width = 1u << bcf_type_shift[val_type_];

      local_data_.resize(width * size_);
      val_ptr_ = local_data_.data();

      switch (val_type_)
      {
      case 0x01u:
        std::copy_n(vec.begin(), size_, (std::int8_t*)val_ptr_);
        break;
      case 0x02u:
        std::copy_n(vec.begin(), size_, (std::int16_t*)val_ptr_); // TODO: handle endianess
        break;
      case 0x03u:
        std::copy_n(vec.begin(), size_, (std::int32_t*)val_ptr_);
        break;
      case 0x04u:
        std::copy_n(vec.begin(), size_, (std::int64_t*)val_ptr_);
        break;
      case 0x05u:
        std::copy_n(vec.begin(), size_, (float*)val_ptr_);
        break;
      }
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, ::savvy::compressed_vector<typename T::value_type>>::value && std::is_signed<typename T::value_type>::value, void>::type
    typed_value::init(const T& vec)
    {
      std::size_t offset_max = 0;
      std::size_t last_off = 0;
      for (auto it = vec.begin(); it != vec.end(); ++it)
      {
        std::size_t off = it.offset() - last_off;
        last_off = it.offset() + 1;
        if (off > offset_max)
          offset_max = off;
      }

      off_type_ = type_code_ignore_missing(static_cast<std::int64_t>(offset_max));

      typedef typename T::value_type vtype;
      if (std::is_integral<vtype>::value && !std::is_same<std::int8_t, vtype>::value)
      {
        vtype min_val = 0;
        vtype max_val = 0;
        for (auto it = vec.begin(); it != vec.end(); ++it)
        {
          if (!is_missing(*it))
          {
            if (*it > max_val)
              max_val = *it;
            if (*it < min_val)
              min_val = *it;
          }
        }

        val_type_ = type_code(std::min(vtype(-max_val), min_val));
      }
      else
      {
        val_type_ = type_code<vtype>();
      }

      sparse_size_ = vec.non_zero_size();
      size_ = vec.size();
      std::size_t off_width = 1u << bcf_type_shift[off_type_];
      std::size_t val_width = 1u << bcf_type_shift[val_type_];
      std::size_t pair_width = off_width + val_width;

      local_data_.resize(pair_width * sparse_size_);
      off_ptr_ = local_data_.data();
      val_ptr_ = local_data_.data() + (sparse_size_ * off_width);

      switch (off_type_)
      {
      case 0x01u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int8_t*)off_ptr_);
        break;
      case 0x02u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int16_t*)off_ptr_); // TODO: handle endianess
        break;
      case 0x03u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int32_t*)off_ptr_);
        break;
      case 0x04u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int64_t*)off_ptr_);
        break;
      }

      switch (val_type_)
      {
      case 0x01u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int8_t*)val_ptr_);
        break;
      case 0x02u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int16_t*)val_ptr_); // TODO: handle endianess
        break;
      case 0x03u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int32_t*)val_ptr_);
        break;
      case 0x04u:
        std::copy_n(vec.value_data(), sparse_size_, (std::int64_t*)val_ptr_);
        break;
      case 0x05u:
        std::copy_n(vec.value_data(), sparse_size_, (float*)val_ptr_);
        break;
      }
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, std::string>::value, void>::type
    typed_value::init(const T& vec)
    {
      val_type_ = typed_value::str;

      size_ = vec.size();

      local_data_.resize(size_);
      val_ptr_ = local_data_.data();
      std::copy_n(vec.begin(), size_, local_data_.begin());
    }

    inline
    site_info::site_info(std::string chrom,
      std::uint32_t pos,
      std::string ref,
      std::vector<std::string> alts,
      std::string id,
      float qual,
      std::vector<std::string> filters,
      std::vector<std::pair<std::string, typed_value>> info)
      :
      chrom_(std::move(chrom)),
      id_(std::move(id)),
      pos_(pos),
      qual_(qual),
      ref_(std::move(ref)),
      alts_(std::move(alts)),
      filters_(std::move(filters)),
      info_(std::move(info))
    {

    }

    inline
    bool site_info::internal::read(site_info& s, std::istream& ifs, const dictionary& dict, std::uint32_t shared_sz)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample

      if (ifs.read((char*)buf.data(), buf.size() * 4))
      {
        std::int32_t tmp_int = buf[0].i;

        if (dict.int_to_str[dictionary::contig].size() <= tmp_int)
        {
          std::fprintf(stderr, "Error: Invalid contig id");
          return false;
        }
        s.chrom_ = dict.int_to_str[dictionary::contig][tmp_int];

        s.pos_ = static_cast<std::uint32_t>(buf[1].i);
        // skip rlen
        s.qual_ = buf[3].f;

        std::uint32_t tmp_uint = static_cast<std::uint32_t>(buf[4].i);
        std::size_t n_allele = tmp_uint >> 16u;
        std::size_t n_info = 0xFFFFu & tmp_uint;

        tmp_uint = static_cast<std::uint32_t>(buf[5].i);
        s.n_fmt_ = tmp_uint >> 24u;
        //std::size_t n_sample = 0xFFFFFFu & tmp_uint;

        s.shared_data_.resize(shared_sz - 6 * 4);
        ifs.read(s.shared_data_.data(), s.shared_data_.size());

        auto shared_it = s.shared_data_.begin();

        try
        {
          // Parse ID
          shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), s.id_);

          // Parse REF/ALT
          if (n_allele)
          {
            shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), s.ref_);
            s.alts_.resize(n_allele - 1);
            for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
            {
              shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), *it);
            }
          }

          // Parse FILTER
          std::vector<std::int32_t> filter_ints;
          shared_it = bcf::deserialize_vec(shared_it, s.shared_data_.end(), filter_ints);
          s.filters_.clear();
          s.filters_.reserve(filter_ints.size());
          for (auto it = filter_ints.begin(); it != filter_ints.end(); ++it)
          {
            if (dict.int_to_str[dictionary::id].size() <= *it)
            {
              std::fprintf(stderr, "Error: Invalid filter id");
              return false;
            }
            s.filters_.emplace_back(dict.int_to_str[dictionary::id][*it]);
          }

          // Parse INFO
          s.info_.resize(n_info);
          auto info_it = s.info_.begin();
          for ( ; info_it != s.info_.end(); ++info_it)
          {
            std::int32_t info_key_id;
            shared_it = bcf::deserialize_int(shared_it, s.shared_data_.end(), info_key_id);
            if (dict.int_to_str[dictionary::id].size() <= info_key_id)
            {
              std::fprintf(stderr, "Error: Invalid info id");
              return false;
            }
            std::string info_key = dict.int_to_str[dictionary::id][info_key_id];

            if (shared_it == s.shared_data_.end())
              break;

            // ------------------------------------------- //
            // TODO: potentially move this to static method since it's similar to FMT parsing
            std::uint8_t type_byte = *(shared_it++);
            std::size_t type_width = 1u << bcf_type_shift[type_byte & 0x0Fu];
            std::size_t sz = type_byte >> 4u;
            if (sz == 15u)
              shared_it = bcf::deserialize_int(shared_it, s.shared_data_.end(), sz);

            if (s.shared_data_.end() - shared_it < (sz * type_width))
              break;

            *info_it = std::make_pair(std::move(info_key), typed_value(type_byte & 0x0Fu, sz, sz * type_width ? &(*shared_it) : nullptr));
            shared_it += sz * type_width;
            // ------------------------------------------- //

          }

          if (info_it == s.info_.end())
            return true;
        }
        catch (const std::exception& e)
        {
          std::fprintf(stderr, "Error: Invalid record data");
          return false;
        }
      }

      std::fprintf(stderr, "Error: Invalid record data");
      return false;
    }

    template <typename Itr>
    bool site_info::internal::serialize(const site_info& s, Itr out_it, const dictionary& dict, std::uint32_t n_sample, std::uint16_t n_fmt)
    {
      union u
      {
        std::int32_t i;
        float f;
      };

      std::array<u, 6> buf; // chrom through n.fmt.sample
      auto res = dict.str_to_int[dictionary::contig].find(s.chrom_);
      if (res == dict.str_to_int[dictionary::contig].end())
      {
        std::fprintf(stderr, "Error: Contig not in header");
        return false;
      }

      buf[0].i = static_cast<std::int32_t>(res->second);
      buf[1].i = static_cast<std::int32_t>(s.pos_);
      buf[2].i = static_cast<std::int32_t>(s.chrom_.size());
      buf[3].f = s.qual_;

      std::uint32_t tmp_uint = (std::uint32_t(s.alts_.size() + 1) << 16u) | (0xFFFFu & std::uint32_t(s.info_.size()));
      buf[4].i = static_cast<std::int32_t>(tmp_uint);

      tmp_uint = (n_fmt << 24u) | (0xFFFFFFu & n_sample);
      buf[5].i = static_cast<std::int32_t>(tmp_uint);

      std::copy_n((char*)buf.data(), buf.size() * sizeof(u), out_it);

      // Encode REF/ALTS
      bcf::serialize_typed_str(out_it, s.id_);
      bcf::serialize_typed_str(out_it, s.ref_);
      for (auto it = s.alts_.begin(); it != s.alts_.end(); ++it)
        bcf::serialize_typed_str(out_it, *it);

      // Encode FILTER
      std::vector<std::int32_t> filter_ints;
      filter_ints.reserve(s.filters_.size());
      for (auto it = s.filters_.begin(); it != s.filters_.end(); ++it)
      {
        res = dict.str_to_int[dictionary::id].find(*it);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: Filter not in header");
          return false;
        }
        filter_ints.emplace_back(res->second);
      }
      bcf::serialize_typed_vec(out_it, filter_ints);

      // Encode INFO
      for (auto it = s.info_.begin(); it != s.info_.end(); ++it)
      {
        res = dict.str_to_int[dictionary::id].find(it->first);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: INFO key not in header");
          return false;
        }

        bcf::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));
        typed_value::internal::serialize(it->second, out_it);
      }

      return true;
    }

    inline
    bool variant::internal::read(variant& v, std::istream& ifs, const dictionary& dict, std::size_t sample_size, bool is_bcf)
    {
      std::uint32_t shared_sz, indiv_sz;
      ifs.read((char*)&shared_sz, sizeof(shared_sz)); // TODO: endian
      ifs.read((char*)&indiv_sz, sizeof(indiv_sz));

      if (ifs && site_info::internal::read(v, ifs, dict, shared_sz))
      {
        v.indiv_buf_.resize(indiv_sz);
        if (ifs.read(v.indiv_buf_.data(), v.indiv_buf_.size()))
        {
          auto indiv_it = v.indiv_buf_.begin();
          v.format_fields_.resize(v.n_fmt_);

          auto fmt_it = v.format_fields_.begin();
          for (; fmt_it != v.format_fields_.end(); ++fmt_it)
          {
            try
            {
              std::int32_t fmt_key_id;
              indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), fmt_key_id);
              if (dict.int_to_str[dictionary::id].size() <= fmt_key_id)
              {
                std::fprintf(stderr, "Error: Invalid FMT id");
                return false;
              }
              std::string fmt_key = dict.int_to_str[dictionary::id][fmt_key_id];

              if (indiv_it == v.indiv_buf_.end())
                break;

              // ------------------------------------------- //
              // TODO: potentially move this to static method since it's similar to INFO parsing.
              std::uint8_t type_byte = *(indiv_it++);
              std::uint8_t type = 0x0Fu & type_byte;

              std::size_t sz = type_byte >> 4u; // TODO: support BCF vector size.
              if (sz == 15u)
                indiv_it = bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), sz);

              if (type == typed_value::sparse)
              {
                assert(!is_bcf);

                if (indiv_it == v.indiv_buf_.end())
                  break;

                std::uint8_t sp_type_byte = *(indiv_it++);
                std::uint8_t off_type = sp_type_byte >> 4u;
                std::uint8_t val_type = sp_type_byte & 0x0Fu;
                std::size_t sp_sz;
                bcf::deserialize_int(indiv_it, v.indiv_buf_.end(), sp_sz);
                std::size_t pair_width = 1u << bcf_type_shift[off_type];
                pair_width += 1u << bcf_type_shift[val_type];

                if (v.indiv_buf_.end() - indiv_it < (sp_sz * pair_width))
                  break;

                *fmt_it = std::make_pair(std::move(fmt_key), typed_value(val_type, sz, off_type, sp_sz, sp_sz * pair_width ? &(*indiv_it) : nullptr));
                indiv_it += sp_sz * pair_width;
              }
              else
              {
                if (is_bcf)
                  sz = sz * sample_size;

                std::size_t type_width = 1u << bcf_type_shift[type];

                if (v.indiv_buf_.end() - indiv_it < (sz * type_width))
                  break;

                *fmt_it = std::make_pair(std::move(fmt_key), typed_value(type, sz, sz * type_width ? &(*indiv_it) : nullptr));
                indiv_it += sz * type_width;
                // ------------------------------------------- //
              }
            }
            catch (const std::exception& e)
            {
              std::fprintf(stderr, "Error: Invalid record data");
              return false;
            }
          }

          if (fmt_it == v.format_fields_.end())
            return true;
        }
      }

      std::fprintf(stderr, "Error: Invalid record data");
      return false;
    }

    inline
    bool variant::internal::serialize(const variant& v, std::vector<char>& buf, const dictionary& dict, std::size_t sample_size, bool is_bcf, std::uint32_t& shared_sz, std::uint32_t& indiv_sz)
    {
      buf.clear();
      buf.reserve(24);

      if (!site_info::internal::serialize(v, std::back_inserter(buf), dict, sample_size, v.format_fields_.size()))
        return false;

      if (buf.size() > std::numeric_limits<std::uint32_t>::max())
      {
        fprintf(stderr, "Error: shared data too big");
        return false;
      }

      shared_sz = buf.size();

      // Encode FMT
      auto out_it = std::back_inserter(buf);
      for (auto it = v.format_fields_.begin(); it != v.format_fields_.end(); ++it)
      {
        auto res = dict.str_to_int[dictionary::id].find(it->first);
        if (res == dict.str_to_int[dictionary::id].end())
        {
          std::fprintf(stderr, "Error: FMT key not in header");
          return false;
        }

        // TODO: Allow for BCF writing.
        bcf::serialize_typed_scalar(out_it, static_cast<std::int32_t>(res->second));
        typed_value::internal::serialize(it->second, out_it);
      }

      if (buf.size() - shared_sz > std::numeric_limits<std::uint32_t>::max())
      {
        std::fprintf(stderr, "Error: individual data too big");
        return false;
      }

      indiv_sz = buf.size() - shared_sz;

      return true;
    }

    template <typename T>
    void variant::set_format(const std::string& key, const std::vector<T>& geno, std::set<std::string> sparse_keys)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          it->second = typed_value(geno);
          return;
        }
      }

      if (it == format_fields_.end())
      {
        format_fields_.emplace_back(key, geno);
      }

      //serialize(format_fields_[it - format_fields_.begin()], geno, sparse_keys.find(key) != sparse_keys.end());
    }

    template <typename T>
    void variant::set_format(const std::string& key, const compressed_vector<T>& geno)
    {
      auto it = format_fields_.begin();
      for ( ; it != format_fields_.end(); ++it)
      {
        if (it->first == key)
        {
          it->second = typed_value(geno);
          return;
        }
      }

      if (it == format_fields_.end())
      {
        format_fields_.emplace_back(key, geno);
      }
    }

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