
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
        this->property_fields_ = {"QUAL", "FILTER"};
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

    reader::reader(const std::string& file_path) :
      hts_file_(bcf_open(file_path.c_str(), "r")),
      hts_hdr_(nullptr),
      hts_rec_(bcf_init1())
    {
      if (!hts_file_ || !hts_rec_)
      {
        this->state_ = std::ios::badbit;
      }
      else
      {
        hts_hdr_ = bcf_hdr_read(hts_file_);
        this->init_property_fields();

        for (int i = 0; i < hts_hdr_->n[BCF_DT_CTG]; ++i)
        {
          auto a = hts_hdr_->id[BCF_DT_CTG][i].key;
          auto b = a;
        }
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

    std::vector<std::string> reader::chromosomes() const
    {
      std::vector<std::string> ret(hts_hdr_->n[BCF_DT_CTG] > 0 ? (unsigned)hts_hdr_->n[BCF_DT_CTG] : 0);
      for (int i = 0; i < ret.size(); ++i)
      {
        ret[i] = hts_hdr()->id[BCF_DT_CTG][i].key;
      }
      return ret;
    }


    bool reader::read_hts_record()
    {
      if (bcf_read(hts_file_, hts_hdr_, hts_rec_) >= 0)
      {
        return true;
      }
      return false;
    }

    indexed_reader::indexed_reader(const std::string& file_path, const region& reg) :
      file_path_(file_path),
      synced_readers_(bcf_sr_init()),
      hts_rec_(nullptr)
    {
      std::stringstream contigs;
      contigs << reg.chromosome();
      if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
        contigs << ":" << reg.from() << "-" << reg.to();

      bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0);
      bcf_sr_add_reader(synced_readers_, file_path_.c_str());
      this->init_property_fields();
    }

    indexed_reader::~indexed_reader()
    {
      if (synced_readers_)
        bcf_sr_destroy(synced_readers_);
    }

    void indexed_reader::reset_region(const region& reg)
    {
      if (synced_readers_)
        bcf_sr_destroy(synced_readers_);
      synced_readers_ = bcf_sr_init();
      state_ = std::ios::goodbit;

      std::stringstream contigs;
      contigs << reg.chromosome();
      if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
        contigs << ":" << reg.from() << "-" << reg.to();

      bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0);
      bcf_sr_add_reader(synced_readers_, file_path_.c_str());
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
