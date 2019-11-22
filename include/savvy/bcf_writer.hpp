#ifndef LIBSAVVY_BCF_WRITER_HPP
#define LIBSAVVY_BCF_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <set>
#include <algorithm>

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
    struct typed_value
    {
      static const std::uint8_t int8   = 1;
      static const std::uint8_t int16  = 2;
      static const std::uint8_t int32  = 3;
      static const std::uint8_t int64  = 4;
      static const std::uint8_t real   = 5;
      static const std::uint8_t real64 = 6;
      static const std::uint8_t str    = 7;
      static const std::uint8_t sparse = 8;

      std::uint8_t type;
      std::size_t length;
      std::vector<std::uint8_t> data;
    };

    class site_info
    {
    private:
      std::string chrom_;
      std::string id_;
      std::int32_t pos_;
      float qual_;
      std::string ref_;
      std::vector<std::string> alts_;
      std::vector<std::string> filters_;
      std::vector<std::pair<std::string, typed_value>> info_;
    public:
      site_info(std::string chrom,
        std::string id,
        std::int32_t pos,
        float qual,
        std::string ref,
        std::vector<std::string> alts,
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

      const std::string& chrom() const { return chrom_; }
      const std::string& id() const { return id_; }
      std::int32_t pos() const { return pos_; }
      float qual() const { return qual_; }
      const std::string& ref() const { return ref_; }
      const std::vector<std::string>& alts() const { return alts_; }
      const std::vector<std::string>& filters() const { return filters_; }
      const std::vector<std::pair<std::string, typed_value>>& info() const { return info_; }
    };

    class variant : public site_info
    {
    public:
      class individual_data
      {
      private:
        std::string key_;
        std::vector<std::uint8_t> data_;
      public:
        individual_data(std::string k, std::vector<std::uint8_t> d) :
          key_(std::move(k)),
          data_(std::move(d))
        {}

        const std::string& key() const { return key_; }
        const std::vector<std::uint8_t>& serialized_data() const { return data_; }
      };
    private:
      std::vector<individual_data> indiv_data_;
    public:
      const std::vector<individual_data>& format_fields() const { return indiv_data_; }
      template <typename T>
      void set_format(const std::string& key, const std::vector<T>& geno, std::set<std::string> sparse_keys = {"GT","EC","DS","HDS"})
      {
        auto it = indiv_data_.begin();
        for ( ; it != indiv_data_.end(); ++it)
        {
          if (it->key() == key)
            break;
        }

        if (it == indiv_data_.end())
        {
          indiv_data_.emplace_back(key, std::vector<std::uint8_t>());
        }

        serialize(indiv_data_[std::distance(indiv_data_.begin(), it)], geno, sparse_keys.find(key) != sparse_keys.end());
      }

      template <typename T>
      void set_format(const std::string& key, const compressed_vector<T>& geno)
      {
        auto it = indiv_data_.begin();
        for ( ; it != indiv_data_.end(); ++it)
        {
          if (it->key() == key)
            break;
        }

        if (it == indiv_data_.end())
        {
          indiv_data_.emplace_back(key, std::vector<std::uint8_t>());
        }

        serialize(indiv_data_[std::distance(indiv_data_.begin(), it)], geno);
      }

    private:
      template <typename T>
      static void serialize_format(std::vector<std::uint8_t>& dest, const std::vector<T>& geno, bool sparse)
      {

      }

      template <typename T>
      static void serialize_format(/* TODO: make this a back_inserter */ std::vector<std::uint8_t>& dest, const compressed_vector<T>& geno)
      {
        std::size_t offset_max = 0;
        std::size_t last_off = 0;
        for (auto it = geno.begin(); it != geno.end(); ++it)
        {
          std::size_t off = it.offset() - last_off;
          last_off = it.offset() + 1;
          if (off > offset_max)
            offset_max = off;
        }

        dest.clear();
        dest.reserve(1 + 8 + geno.non_zero_size() + sizeof(T) * geno.non_zero_size()); // TODO make this more accurate.

        std::uint8_t type_byte = typed_value::sparse;
        if (geno.size() >= 15u)
          type_byte = (15u << 4u) | type_byte;
        else
          type_byte = (geno.size() << 4u) | type_byte;
        dest.emplace_back(type_byte);

        if (geno.size() >= 15u)
        {
          bcf::write_typed_scalar(std::back_inserter(dest), geno.non_zero_size());
        }

        // TODO: write off_type|val_type, sparse size, and data.
      }
    };


//    struct packed_record
//    {
//      std::uint32_t l_shared = 0;
//      std::uint32_t l_indiv = 0;
//      std::int32_t chrom_idx = 1;
//      std::int32_t pos = 123;
//      std::int32_t ref_len = 1;
//      float qual = 1.f;
//      std::uint32_t n_allele_info = (2 << 16) | 0;
//      std::uint32_t n_fmt_sample = (1 << 24) | 6;
//      std::string variant_id = "rs123";
//      std::vector<std::string> alleles = {"C","T"};
//      std::vector<std::string> filters;
//      std::vector<std::pair<std::string, typed_value>> info;
//    };

    class dictionary
    {
    public:
      static const std::uint8_t id = 0;
      static const std::uint8_t contig = 1;
      static const std::uint8_t sample = 2;
      std::array<std::unordered_map<std::string, std::uint32_t>, 3> str_to_int;
      std::array<std::vector<std::string>, 3> int_to_str;
    };


    void write_header(std::ostream& os, const std::vector<std::pair<std::string, std::string>>& headers, const std::vector<std::string>& ids)
    {
      std::string magic = "SAV\x02\x00";

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

      os.write(magic.data(), magic.size());
      os.write((char*)(&header_block_sz), sizeof(header_block_sz));

      for (auto it = headers.begin(); it != headers.end(); ++it)
      {
        os.write("##", 2);
        os.write(it->first.data(), it->first.size());
        os.write("=", 1);
        os.write(it->second.data(), it->second.size());
        os.write("\n", 1);
      }

      os.write(column_names.data(), column_names.size());
      for (auto it = ids.begin(); it != ids.end(); ++it)
      {
        os.write("\t", 1);
        os.write(it->data(), it->size());
      }

      os.write("\n\0", 2);
    }

    void write_record(std::ostream& os, const dictionary& dict, const variant& r)
    {
      std::uint32_t l_shared = 6 * 4; // chrom through n.fmt.sample

      l_shared += bcf::get_typed_value_size(r.id());

      l_shared += bcf::get_typed_value_size(r.ref());
      for (auto& a : r.alts())
        l_shared += bcf::get_typed_value_size(a);

      std::vector<std::int8_t> filter_vec(r.filters().size()); // TODO: Allow more than 255 filters.
      for (std::size_t i = 0; i < filter_vec.size(); ++i)
      {
        auto it = dict.str_to_int[dictionary::id].find(r.filters()[i]);
        if (it == dict.str_to_int[dictionary::id].end())
        {
          std::cerr << "Error: filter not valid" << std::endl;
          os.setstate(os.rdstate() | std::ios::badbit);
        }
        assert(it->second <= 0xFF);
        filter_vec[i] = it->second;
      }

      l_shared += bcf::get_typed_value_size(filter_vec);
      for (auto& i: r.info())
      {
        auto it = dict.str_to_int[dictionary::contig].find(i.first);
        if (it == dict.str_to_int[dictionary::contig].end())
        {
          std::cerr << "Error: INFO key not valid" << std::endl;
          os.setstate(os.rdstate() | std::ios::badbit);
          return;
        }

        l_shared += bcf::get_typed_value_size((std::int32_t)it->second); // TODO: Allow for variable sized INFO keys.
        l_shared += i.second.data.size(); //(bcf::get_typed_value_size(0) + bcf::get_typed_value_size(i.second));
      }


      std::uint32_t l_indiv = 0;
      for (auto& f : r.format_fields())
      {
        auto it = dict.str_to_int[dictionary::id].find(f.key());
        if (it != dict.str_to_int[dictionary::id].end())
        {
          l_indiv += bcf::get_typed_value_size((std::int32_t)it->second);
          l_indiv += f.serialized_data().size();
        }
      }

      os.write((char*)(&l_shared), sizeof(l_shared));
      os.write((char*)(&l_indiv), sizeof(l_indiv));

      {
        auto it = dict.str_to_int[dictionary::contig].find(r.chrom());
        if (it == dict.str_to_int[dictionary::contig].end())
        {
          std::cerr << "Error: contig not valid" << std::endl;
          os.setstate(os.rdstate() | std::ios::badbit);
          return;
        }
        os.write((char*)(&it->second), sizeof(it->second));
      }


      std::uint32_t tmp_i32 = r.pos();
      os.write((char*)(&tmp_i32), sizeof(tmp_i32));
      tmp_i32 = r.ref().size();
      os.write((char*)(&tmp_i32), sizeof(tmp_i32));
      float qual = r.qual();
      os.write((char*)(&qual), sizeof(qual));

      tmp_i32 = ((1 + r.alts().size()) << 16) | r.info().size();
      os.write((char*)(&tmp_i32), sizeof(tmp_i32));
      tmp_i32 = (r.format_fields().size() << 24) | 0; // zero for n_sample in SAV
      os.write((char*)(&tmp_i32), sizeof(tmp_i32));


      bcf::write_typed_str(os, r.id());
      bcf::write_typed_str(os, r.ref());
      for (auto& a: r.alts())
        bcf::write_typed_str(os, a);
      bcf::write_typed_vec(os, filter_vec);

      for (auto& i: r.info())
      {
        bcf::write_typed_scalar(os, (std::int32_t)dict.str_to_int[dictionary::id].find(i.first)->second); // TODO: Allow for variable sized INFO keys.
        bcf::write_type_and_size(os, i.second.type, i.second.length);
        os.write((char*)(i.second.data.data()), i.second.data.size());
      }

      for (auto& f : r.format_fields())
      {
        bcf::write_typed_scalar(os, (std::int32_t)(dict.str_to_int[dictionary::id].find(f.key())->second)); // TODO: Allow for variable sized INFO keys.
        os.write((char*)(f.serialized_data().data()), f.serialized_data().size());
      }
    }
  }
}

#endif // LIBSAVVY_BCF_WRITER_HPP