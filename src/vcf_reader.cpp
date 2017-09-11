
#include "savvy/vcf_reader.hpp"


#include <htslib/vcf.h>

#include <sstream>
#include <assert.h>
#include <limits>

namespace savvy
{
  namespace vcf
  {

    const char** reader_base::samples_begin() const
    {
      return hts_hdr() ? (const char**)hts_hdr()->samples : nullptr;
    }

    const char** reader_base::samples_end() const
    {
      return hts_hdr() ? (const char**)((hts_hdr()->samples) + bcf_hdr_nsamples(hts_hdr())) : nullptr;
    }

    std::uint64_t reader_base::sample_count() const
    {
      return static_cast<std::uint64_t>(bcf_hdr_nsamples(hts_hdr()));
    }

    std::vector<std::pair<std::string, std::string>> reader_base::headers() const
    {
      std::vector<std::pair<std::string, std::string>> ret;

      bcf_hdr_t* hdr = hts_hdr();
      if (hdr)
      {
        ret.reserve(hdr->nhrec - 1);
        for (int i = 1; i < hdr->nhrec; ++i)
        {
          std::string key, val;
          if (hdr->hrec[i]->key && hdr->hrec[i]->value)
          {
            key = hdr->hrec[i]->key;
            val = hdr->hrec[i]->value;
          }
          else if (hdr->hrec[i]->key &&  (hdr->hrec[i]->type == BCF_HL_INFO || hdr->hrec[i]->type == BCF_HL_FLT || hdr->hrec[i]->type == BCF_HL_STR))
          {
            bcf_hrec_t* r = hdr->hrec[i];
            key = r->key;
            std::stringstream ss_val;

            ss_val << "<";
            for (int j = 0; j < r->nkeys - 1; ++j) // minus 1 to remove IDX;
            {
              if (j > 0)
                ss_val << ",";
              if (r->keys[j])
                ss_val << r->keys[j];
              ss_val << "=";
              if (r->vals[j])
                ss_val << r->vals[j];
            }
            ss_val << ">";
            val = ss_val.str();
          }

          if (key.size())
            ret.emplace_back(std::move(key), std::move(val));
            //ret.insert(std::upper_bound(ret.begin(), ret.end(), std::make_pair(key, std::string()), [](const auto& a, const auto& b) { return a.first < b.first; }), {std::move(key), std::move(val)});
        }
      }

      return ret;
    }

    std::vector<std::string> reader_base::prop_fields() const
    {
      std::vector<std::string> ret(property_fields_);
      return ret;
    }

    void reader_base::init_property_fields()
    {
      bcf_hdr_t* hdr = hts_hdr();
      if (hdr)
      {
        this->property_fields_ = {"ID","QUAL", "FILTER"};
        for (int i = 0; i < hdr->nhrec; ++i)
        {
          if (hdr->hrec[i]->type == BCF_HL_INFO)
          {
            bcf_hrec_t* r = hdr->hrec[i];
            for (int j = 0; j < r->nkeys; ++j)
            {
              if (strcmp(r->keys[j], "ID") == 0)
              {
                const char* inf = r->vals[j];
                if (inf)
                  this->property_fields_.emplace_back(inf);
              }
            }
          }
        }
      }
    }

    reader::reader(const std::string& file_path, fmt data_format) :
      hts_file_(bcf_open(file_path.c_str(), "r")),
      hts_hdr_(nullptr),
      hts_rec_(bcf_init1())
    {

      requested_data_formats_ = {data_format};
      if (!hts_file_ || !hts_rec_)
      {
        this->state_ = std::ios::badbit;
      }
      else
      {
        hts_hdr_ = bcf_hdr_read(hts_file_);
        this->init_property_fields();
      }
    }

    reader::reader(reader&& source) :
      reader_base(std::move(source)),
      hts_file_(source.hts_file_),
      hts_hdr_(source.hts_hdr_),
      hts_rec_(source.hts_rec_)
    {
      source.hts_file_ = nullptr;
      source.hts_hdr_ = nullptr;
      source.hts_rec_ = nullptr;
      source.gt_ = nullptr;
      source.state_ = std::ios::badbit;
    }

    reader::~reader()
    {
      try
      {
        if (hts_hdr_)
          bcf_hdr_destroy(hts_hdr_);
        if (hts_file_)
          bcf_close(hts_file_);
        if (hts_rec_)
          bcf_destroy1(hts_rec_);
      }
      catch (...)
      {
        assert(!"bcf_hdr_destroy or bcf_close threw exception!");
      }
    }

//    std::vector<std::string> reader::chromosomes() const
//    {
//      std::vector<std::string> ret(hts_hdr_->n[BCF_DT_CTG] > 0 ? (unsigned)hts_hdr_->n[BCF_DT_CTG] : 0);
//      for (int i = 0; i < ret.size(); ++i)
//      {
//        ret[i] = hts_hdr()->id[BCF_DT_CTG][i].key;
//      }
//      return ret;
//    }


    bool reader::read_hts_record()
    {
      if (bcf_read(hts_file_, hts_hdr_, hts_rec_) >= 0)
      {
        return true;
      }
      return false;
    }

    indexed_reader::indexed_reader(const std::string& file_path, const region& reg, fmt data_format) :
      file_path_(file_path),
      synced_readers_(bcf_sr_init()),
      hts_rec_(nullptr)
    {
      this->requested_data_formats_ = {data_format};
      std::stringstream contigs;
      contigs << reg.chromosome();
      if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
        contigs << ":" << reg.from() << "-" << reg.to();

      bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0);
      if (bcf_sr_add_reader(synced_readers_, file_path_.c_str()))
        this->init_property_fields();
      else
        state_ = std::ios::badbit;
    }

    indexed_reader::~indexed_reader()
    {
      if (synced_readers_)
        bcf_sr_destroy(synced_readers_);
    }

    std::vector<std::string> indexed_reader::chromosomes() const
    {
      std::vector<std::string> ret;

      if (good())
      {
        hts_idx_t* idx = bcf_index_load(file_path_.c_str());
        if (idx)
        {
          int n{};
          const char** arr = bcf_index_seqnames(idx, hts_hdr(), &n);
          ret.resize(n);
          for (int i = 0; i < n; ++i)
          {
            ret[i] = arr[i];
          }
        }
      }

      return ret;
    }

    void indexed_reader::reset_region(const region& reg)
    {
      if (good())
      {
        if (synced_readers_)
          bcf_sr_destroy(synced_readers_);
        synced_readers_ = bcf_sr_init();
        hts_rec_ = nullptr;
        state_ = std::ios::goodbit;

        std::stringstream contigs;
        contigs << reg.chromosome();
        if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
          contigs << ":" << reg.from() << "-" << reg.to();

        if (bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0) != 0 || bcf_sr_add_reader(synced_readers_, file_path_.c_str()) != 1)
          state_ = std::ios::failbit;
      }
    }

    bool indexed_reader::read_hts_record()
    {
      if (bcf_sr_next_line(synced_readers_) && (hts_rec_ = bcf_sr_get_line(synced_readers_, 0)))
      {
        return true;
      }
      return false;
    }
  }
}

//savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a)
//{
//  savvy::vcf::marker::const_iterator ret(a);
//  return (ret += n);
//}
