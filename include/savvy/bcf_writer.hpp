#ifndef LIBSAVVY_BCF_WRITER_HPP
#define LIBSAVVY_BCF_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>

namespace savvy
{
  namespace bcf
  {

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint32_t>::type
    get_typed_value_size(std::vector<T> vec)
    {
      static_assert(!std::is_same<T, std::int64_t>::value, "64-bit integers not allowed in BCF spec.");
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

        ret += vec.size() * sizeof(T);
      return ret;
    }

    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, std::uint32_t>::type
    get_typed_value_size(T vec)
    {
      static_assert(!std::is_same<T, std::int64_t>::value, "64-bit integers not allowed in BCF spec.");
      return 1 + sizeof(T);
    }


    struct record
    {
      std::uint32_t l_shared;
      std::uint32_t l_indiv;
      std::int32_t chrom_idx;
      std::int32_t pos;
      std::int32_t ref_len;
      float qual;
      std::uint32_t n_allele_info;
      std::uint32_t n_fmt_sample;

    };

    inline void write_bcf()
    {
      std::ofstream ofs("test.bcf", std::ios::binary);

      // BEGIN HEADER
      std::string magic = "BCF\2\1";

      std::vector<std::pair<std::string, std::string>> headers;
      std::uint32_t header_block_sz = 0;
      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        header_block_sz += it->first.size();
        header_block_sz += it->second.size();
        header_block_sz += 4;
      }

      std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tFILTER\tINFO\tFORMAT";
      header_block_sz += column_names.size();

      std::vector<std::string> samples;
      for (auto it = samples.begin(); it != samples.end(); ++it)
      {
        header_block_sz += 1 + it->size();
      }

      header_block_sz += 2; //new line and null

      ofs.write(magic.data(), magic.size());
      ofs.write((char*)(&header_block_sz), sizeof(header_block_sz));

      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        ofs.write("##", 2);
        ofs.write(it->first.data(), it->first.size());
        ofs.write("=", 1);
        ofs.write(it->second.data(), it->second.size());
        ofs.write("\n", 1);
      }

      ofs.write(column_names.data(), column_names.size());
      for (auto it = samples.begin(); it != samples.end(); ++it)
      {
        ofs.write("\t", 1);
        ofs.write(it->data(), it->size());
      }

      ofs.write("\n\0", 1);
      // END HEADER

      // BEGIN RECORDS
      std::vector<std::int32_t> test;
      get_typed_value_size(test);
      std::int32_t scalar_test;
      get_typed_value_size(scalar_test);

//      std::int64_t scalar64_test;
//      get_typed_value_size(scalar64_test);
      // END RECORDS

    }
  }
}

#endif // LIBSAVVY_BCF_WRITER_HPP