
#include "vcf_reader.hpp"


#include <sstream>
#include <assert.h>
#include <limits>
#include <dep/htslib_project-prefix/src/htslib_project/htslib/vcf.h>

namespace vc
{
  namespace vcf
  {
//    marker::marker(bcf1_t* hts_rec, int* gt, int num_gt, std::uint16_t allele_index) :
//      hts_rec_(hts_rec),
//      gt_(gt),
//      num_gt_(num_gt),
//      allele_index_(allele_index)
//    {
//
//    }
//
//    marker::~marker()
//    {
//
//    }
//
//    const allele_status marker::const_is_missing = allele_status::is_missing;
//    const allele_status marker::const_has_ref = allele_status::has_ref;
//    const allele_status marker::const_has_alt = allele_status::has_alt;
//
//    std::int32_t marker::chrom_id() const
//    {
//      return hts_rec_->rid;
//    }
//
//    std::uint64_t marker::pos() const
//    {
//      return static_cast<std::uint64_t>(hts_rec_->pos + 1);
//    }
//
//    std::string marker::ref() const
//    {
//      std::string ret(hts_rec_->d.allele[0]);
//      return ret;
//    }
//
//    std::string marker::alt() const
//    {
//      if (hts_rec_->n_allele > 1)
//        return std::string(hts_rec_->d.allele[allele_index_]);
//      return "";
//    }
//
//    const allele_status& marker::operator[](std::size_t i) const
//    {
//      if (gt_[i] == bcf_gt_missing)
//        return const_is_missing;
//      else
//        return (bcf_gt_allele(gt_[i]) == allele_index_ ? const_has_alt : const_has_ref);
//    }
//
//    block::block() :
//      hts_rec_(bcf_init1()),
//      gt_(nullptr),
//      gt_sz_(0),
//      num_samples_(0)
//    {
//
//    }
//
//    block::block(block&& source) :
//      hts_rec_(nullptr),
//      gt_(nullptr),
//      gt_sz_(0),
//      num_samples_(0)
//    {
//      operator=(std::move(source));
//    }
//
//    block& block::operator=(block&& source)
//    {
//      if (hts_rec_)
//        bcf_destroy1(hts_rec_);
//      if (gt_)
//        free(gt_);
//
//      hts_rec_ = source.hts_rec_;
//      source.hts_rec_ = nullptr;
//      gt_ = source.gt_;
//      source.gt_ = nullptr;
//      gt_sz_ = source.gt_sz_;
//      num_samples_ = source.num_samples_;
//      markers_ = std::move(source.markers_);
//
//      return *this;
//    }
//
//    block::~block()
//    {
//      try
//      {
//        if (hts_rec_)
//          bcf_destroy1(hts_rec_);
//        if (gt_)
//          free(gt_);
//      }
//      catch (...)
//      {
//        assert(!"bcf_destroy1 threw exception!");
//      }
//    }
//
//    const marker& block::operator[](std::size_t i) const
//    {
//      return markers_[i];
//    }
//
//    bool block::read_block(block& destination, htsFile* hts_file, bcf_hdr_t* hts_hdr)
//    {
//      destination = block();
//      if (bcf_read(hts_file, hts_hdr, destination.hts_rec_) >= 0)
//      {
//        bcf_unpack(destination.hts_rec_, BCF_UN_ALL); //BCF_UN_STR | BCF_UN_FMT);
//        bcf_get_genotypes(hts_hdr, destination.hts_rec_, &(destination.gt_), &(destination.gt_sz_));
//        destination.num_samples_ = hts_hdr->n[BCF_DT_SAMPLE];
//        if (destination.gt_sz_ % destination.num_samples_ != 0)
//        {
//          // TODO: mixed ploidy at site error.
//        }
//        else
//        {
//          std::uint16_t i = 1;
//          do
//          {
//            destination.markers_.emplace_back(destination.hts_rec_, destination.gt_, destination.gt_sz_, i);
//            ++i;
//          } while (i < destination.hts_rec_->n_allele);
//          return true;
//        }
//      }
//      return false;
//    }
//
//    bool block::read_block(block& destination, bcf_srs_t* sr)
//    {
//      destination = block();
//      bcf1_t* rec;
//      bcf_hdr_t* hdr = bcf_sr_get_header(sr, 0);
//      if (bcf_sr_next_line(sr) && (rec = bcf_sr_get_line(sr, 0)))
//      {
//        bcf_copy(destination.hts_rec_, rec);
//        bcf_unpack(destination.hts_rec_, BCF_UN_ALL); //BCF_UN_STR | BCF_UN_FMT);
//        bcf_get_genotypes(hdr, destination.hts_rec_, &(destination.gt_), &(destination.gt_sz_));
//        destination.num_samples_ = hdr->n[BCF_DT_SAMPLE];
//        if (destination.gt_sz_ % destination.num_samples_ != 0)
//        {
//          // TODO: mixed ploidy at site error.
//        }
//        else
//        {
//          std::uint16_t i = 1;
//          do
//          {
//            destination.markers_.emplace_back(destination.hts_rec_, destination.gt_, destination.gt_sz_, i);
//            ++i;
//          } while (i < destination.hts_rec_->n_allele);
//          return true;
//        }
//      }
//      return false;
//    }

    char** reader_base::samples_begin() const
    {
      return hts_hdr() ? hts_hdr()->samples : nullptr;
    }

    char** reader_base::samples_end() const
    {
      return hts_hdr() ? hts_hdr()->samples + bcf_hdr_nsamples(hts_hdr()) : nullptr;
    }

    std::uint64_t reader_base::sample_count() const
    {
      return static_cast<std::uint64_t>(bcf_hdr_nsamples(hts_hdr()));
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

//    reader& reader::operator>>(block& destination)
//    {
//      if (!block::read_block(destination, hts_file_, hts_hdr_))
//        this->state_ = std::ios::failbit;
//      return *this;
//    }

//    std::string reader::get_chromosome(const reader& rdr, const marker& mkr)
//    {
//      std::string ret(bcf_hdr_id2name(rdr.hts_hdr(), mkr.chrom_id()));
//      return ret;
//    }

//    index_reader::index_reader(const std::string& file_path)
//    {
//      htsFile* f = bcf_open(file_path.c_str(), "r");
//      if (!f)
//      {
//        this->state_ |= std::ios::badbit;
//      }
//      else
//      {
//        bcf_hdr_t* hdr = bcf_hdr_read(f);
//        if (!hdr)
//        {
//          this->state_ |= std::ios::badbit;
//        }
//        else
//        {
//          std::stringstream contigs;
//          for (int i = 0; i < hdr->n[BCF_DT_CTG]; ++i)
//          {
//            if (i > 0)
//              contigs << ",";
//            contigs << hdr->id[BCF_DT_CTG][i].key;
//          }
//
//          synced_readers_ = bcf_sr_init();
//
//          bcf_sr_set_regions(synced_readers_, contigs.str().c_str(), 0);
//          bcf_sr_add_reader(synced_readers_, file_path.c_str());
//          bcf_sr_seek(synced_readers_, "", -1);
//
//
//
//          bcf_hdr_destroy(hdr);
//        }
//        bcf_close(f);
//      }
//    }
//
//    index_reader::~index_reader()
//    {
//
//      if (synced_readers_)
//        bcf_sr_destroy(synced_readers_);
//    }
//
//    index_reader& index_reader::operator>>(block& destination)
//    {
//      if (!block::read_block(destination, synced_readers_))
//        this->state_ = std::ios::failbit;
//      return *this;
//    }
//
//    index_reader& index_reader::seek(const std::string& chromosome, std::uint64_t position)
//    {
//      bcf_sr_seek(synced_readers_, chromosome.c_str(), static_cast<int>(position & std::numeric_limits<int>::max()) - 1);
//      return *this;
//    }
//
//    std::string index_reader::get_chromosome(const index_reader& rdr, const marker& mkr)
//    {
//      std::string ret(bcf_hdr_id2name(rdr.hts_hdr(), mkr.chrom_id()));
//      return ret;
//    }
  }
}

//vc::vcf::marker::const_iterator operator+(vc::vcf::marker::const_iterator::difference_type n, const vc::vcf::marker::const_iterator& a)
//{
//  vc::vcf::marker::const_iterator ret(a);
//  return (ret += n);
//}
