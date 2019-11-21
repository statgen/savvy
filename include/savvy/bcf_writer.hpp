#ifndef LIBSAVVY_BCF_WRITER_HPP
#define LIBSAVVY_BCF_WRITER_HPP

#include "savvy/varint.hpp"

#include <string>
#include <vector>
#include <fstream>

namespace savvy
{
  namespace bcf
  {
    template <typename T>
    typename std::enable_if<std::is_signed<T>::value, void>::type
    write_typed_scalar(std::ostream& os, const T& val)
    {
      std::uint8_t type_byte = 1;
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
      os.write((char*)&val, sizeof(val));
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


    void write_typed_str(std::ostream& os, const std::string& str)
    {
      std::uint8_t type_byte = str.size() < 15 ? str.size() : 15;
      type_byte = (type_byte << 4) | 0x07;


      os.write((char*)&type_byte, 1);
      if (str.size() >= 15)
      {
        if (str.size() <= 0x7F)
          write_typed_scalar(os, (std::int8_t)str.size());
        else if (str.size() <= 0x7FFF)
          write_typed_scalar(os, (std::int16_t)str.size());
        else if (str.size() <= 0x7FFFFFFF)
          write_typed_scalar(os, (std::int32_t)str.size());
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


    struct record
    {
      std::uint32_t l_shared = 0;
      std::uint32_t l_indiv = 0;
      std::int32_t chrom_idx = 1;
      std::int32_t pos = 123;
      std::int32_t ref_len = 1;
      float qual = 1;
      std::uint32_t n_allele_info = (2 << 16) | 0;
      std::uint32_t n_fmt_sample = (1 << 24) | 6;
      std::string variant_id = "rs123";
      std::vector<std::string> alleles = {"C","T"};
      std::vector<std::string> filters;
      std::vector<std::pair<std::string, std::string>> info;
    };

    class writer
    {
    private:
      std::string file_path_;
      std::ofstream file_;
      std::vector<std::pair<std::string, std::string>> headers_;
      std::vector<std::string> samples_;
    public:
      writer(const std::string& file_path,
        std::vector<std::pair<std::string, std::string>> headers,
        std::vector<std::string> samples) :
        file_path_(file_path),
        file_(file_path, std::ios::binary),
        headers_(std::move(headers)),
        samples_(std::move(samples))
      {

      }

      void write_header()
      {
        std::string magic = "BCF\x02\x02";

        std::uint32_t header_block_sz = 0;
        for (auto it = headers_.begin(); it != headers_.end(); ++it)
        {
          header_block_sz += it->first.size();
          header_block_sz += it->second.size();
          header_block_sz += 4;
        }

        std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        header_block_sz += column_names.size();

        std::vector<std::string> samples = {"1","2","3","4","5","6"};
        for (auto it = samples.begin(); it != samples.end(); ++it)
        {
          header_block_sz += 1 + it->size();
        }

        header_block_sz += 2; //new line and null

        file_.write(magic.data(), magic.size());
        file_.write((char*)(&header_block_sz), sizeof(header_block_sz));

        for (auto it = headers_.begin(); it != headers_.end(); ++it)
        {
          file_.write("##", 2);
          file_.write(it->first.data(), it->first.size());
          file_.write("=", 1);
          file_.write(it->second.data(), it->second.size());
          file_.write("\n", 1);
        }

        file_.write(column_names.data(), column_names.size());
        for (auto it = samples.begin(); it != samples.end(); ++it)
        {
          file_.write("\t", 1);
          file_.write(it->data(), it->size());
        }

        file_.write("\n\0", 2);

      }

      writer& write(record r, std::vector<std::int8_t>& genos)
      {
        r.l_shared += sizeof(r.chrom_idx);
        r.l_shared += sizeof(r.pos);
        r.l_shared += sizeof(r.ref_len);
        r.l_shared += sizeof(r.qual);
        r.l_shared += sizeof(r.n_allele_info);
        r.l_shared += sizeof(r.n_fmt_sample);
        r.l_shared += get_typed_value_size(r.variant_id);

        for (auto& a : r.alleles)
          r.l_shared += get_typed_value_size(a);

        std::vector<std::int8_t> filter_vec = {0};
        r.l_shared += get_typed_value_size(filter_vec);
        for (auto& i: r.info)
          r.l_shared += (get_typed_value_size(0) + get_typed_value_size(i.second));

        r.l_indiv += get_typed_value_size(std::int8_t(0));
        r.l_indiv += 1;
        r.l_indiv += genos.size() * sizeof(std::vector<std::int8_t>::value_type);

        file_.write((char*)(&r.l_shared), sizeof(r.l_shared));
        file_.write((char*)(&r.l_indiv), sizeof(r.l_indiv));
        file_.write((char*)(&r.chrom_idx), sizeof(r.chrom_idx));
        file_.write((char*)(&r.pos), sizeof(r.pos));
        file_.write((char*)(&r.ref_len), sizeof(r.ref_len));
        file_.write((char*)(&r.qual), sizeof(r.qual));
        file_.write((char*)(&r.n_allele_info), sizeof(r.n_allele_info));
        file_.write((char*)(&r.n_fmt_sample), sizeof(r.n_fmt_sample));


        write_typed_str(file_, r.variant_id);
        write_typed_str(file_, r.alleles[0]);
        write_typed_str(file_, r.alleles[1]);
        write_typed_vec(file_, filter_vec);
        write_typed_scalar(file_, std::int8_t(3));

        std::uint16_t ploidy = 2;
        std::uint8_t typed_byte = (ploidy << 4u) | 1u;
        file_.write((char*)&typed_byte, 1);
        file_.write((char*)genos.data(), genos.size());

        return *this;
      }
    };

    inline void write_bcf()
    {
      std::ofstream ofs("test.bcf", std::ios::binary);

      // BEGIN HEADER
      std::string magic = "BCF\x02\x02";

      std::vector<std::pair<std::string, std::string>> headers = {
        {"fileformat","VCFv4.3"},
        {"FILTER","<ID=PASS,Description=\"Pass\">"},
        {"FILTER","<ID=q10,Description=\"Quality below 10\">"},
        {"FILTER","<ID=s50,Description=\"Less than 50% of samples have data\">"},
        {"FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">"},
        {"fileDate","20090805"},
        {"source","myImputationProgramV3.1"},
        {"reference","file:///seq/references/1000GenomesPilot-NCBI36.fasta"},
        {"contig","<ID=18,length=80373285,assembly=B36,species=\"Homo sapiens\",taxonomy=x>"},
        {"contig","<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"},
        {"phasing","partial"},
        {"INFO","<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"},
        {"INFO","<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"},
        {"INFO","<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"},
        {"INFO","<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"},
        {"INFO","<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">"},
        {"INFO","<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">"}};


      std::uint32_t header_block_sz = 0;
      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        header_block_sz += it->first.size();
        header_block_sz += it->second.size();
        header_block_sz += 4;
      }

      std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
      header_block_sz += column_names.size();

      std::vector<std::string> samples = {"1","2","3","4","5","6"};
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

      ofs.write("\n\0", 2);
      // END HEADER

      // BEGIN RECORDS
      record r;
      r.l_shared += sizeof(r.chrom_idx);
      r.l_shared += sizeof(r.pos);
      r.l_shared += sizeof(r.ref_len);
      r.l_shared += sizeof(r.qual);
      r.l_shared += sizeof(r.n_allele_info);
      r.l_shared += sizeof(r.n_fmt_sample);
      r.l_shared += get_typed_value_size(r.variant_id);
      for (auto& a : r.alleles)
        r.l_shared += get_typed_value_size(a);
      std::vector<std::int8_t> filter_vec = {0};
      r.l_shared += get_typed_value_size(filter_vec);
      for (auto& i: r.info)
        r.l_shared += (get_typed_value_size(0) + get_typed_value_size(i.second));

      r.l_indiv += get_typed_value_size(std::int8_t(0));
      r.l_indiv += 1;
      r.l_indiv += 6 * 2 * sizeof(char);

      ofs.write((char*)(&r.l_shared), sizeof(r.l_shared));
      ofs.write((char*)(&r.l_indiv), sizeof(r.l_indiv));
      ofs.write((char*)(&r.chrom_idx), sizeof(r.chrom_idx));
      ofs.write((char*)(&r.pos), sizeof(r.pos));
      ofs.write((char*)(&r.ref_len), sizeof(r.ref_len));
      ofs.write((char*)(&r.qual), sizeof(r.qual));
      ofs.write((char*)(&r.n_allele_info), sizeof(r.n_allele_info));
      ofs.write((char*)(&r.n_fmt_sample), sizeof(r.n_fmt_sample));


      write_typed_str(ofs, r.variant_id);
      write_typed_str(ofs, r.alleles[0]);
      write_typed_str(ofs, r.alleles[1]);
      write_typed_vec(ofs, filter_vec);
      write_typed_scalar(ofs, std::int8_t(3));

      std::uint8_t typed_byte = (2 << 4) | 1;
      ofs.write((char*)&typed_byte, 1);
      ofs.write("\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02", 12);




//      std::int64_t scalar64_test;
//      get_typed_value_size(scalar64_test);
      // END RECORDS

    }
  }

  namespace sav2
  {
    void write_header(std::ostream& os, const std::vector<std::pair<std::string, std::string>>& headers, const std::vector<std::string>& ids)
    {
      std::string magic = "SAV\x02\x00";

      std::uint32_t header_block_sz = 0;
      for (auto it = ids.begin(); it != ids.end(); ++it)
      {
        header_block_sz += it->first.size();
        header_block_sz += it->second.size();
        header_block_sz += 4;
      }

      std::string column_names = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
      header_block_sz += column_names.size();

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

      ofs.write("\n\0", 2);
    }

    void write_site_info(std::ostream& os, const record& r)
    {

    }
  }
}

#endif // LIBSAVVY_BCF_WRITER_HPP