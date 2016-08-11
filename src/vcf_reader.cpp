
#include "vcf_reader.hpp"

#include <assert.h>
#include <vcf.h>

namespace vc
{
  namespace vcf
  {
    marker::marker(bcf1_t* hts_rec, int* gt, int num_gt, std::uint16_t allele_index) :
      hts_rec_(hts_rec),
      gt_(gt),
      num_gt_(num_gt),
      allele_index_(allele_index)
    {

    }

    marker::~marker()
    {

    }

    const allele_status marker::const_is_missing = allele_status::is_missing;
    const allele_status marker::const_has_ref = allele_status::has_ref;
    const allele_status marker::const_has_alt = allele_status::has_alt;

    std::int32_t marker::chrom_id() const
    {
      return hts_rec_->rid;
    }

    std::uint64_t marker::pos() const
    {
      return static_cast<std::uint64_t>(hts_rec_->pos + 1);
    }

    std::string marker::ref() const
    {
      std::string ret(hts_rec_->d.allele[0]);
      return ret;
    }

    std::string marker::alt() const
    {
      if (hts_rec_->n_allele > 1)
        return std::string(hts_rec_->d.allele[allele_index_]);
      return "";
    }

    const allele_status& marker::operator[](std::size_t i) const
    {
      if (gt_[i] == bcf_gt_missing)
        return const_is_missing;
      else
        return (bcf_gt_allele(gt_[i]) == allele_index_ ? const_has_alt : const_has_ref);
    }

    block::block() :
      hts_rec_(bcf_init1()),
      gt_(nullptr),
      gt_sz_(0),
      num_samples_(0)
    {

    }

    block::block(block&& source)
    {
      operator=(std::move(source));
    }

    block& block::operator=(block&& source)
    {
      if (hts_rec_)
        bcf_destroy1(hts_rec_);
      if (gt_)
        free(gt_);

      hts_rec_ = source.hts_rec_;
      source.hts_rec_ = nullptr;
      gt_ = source.gt_;
      source.gt_ = nullptr;
      gt_sz_ = source.gt_sz_;
      num_samples_ = source.num_samples_;
      markers_ = std::move(source.markers_);

      return *this;
    }

    block::~block()
    {
      try
      {
        if (hts_rec_)
          bcf_destroy1(hts_rec_);
        if (gt_)
          free(gt_);
      }
      catch (...)
      {
        assert(!"bcf_destroy1 threw exception!");
      }
    }

    const marker& block::operator[](std::size_t i) const
    {
      return markers_[i];
    }

    bool block::read_block(block& destination, htsFile* hts_file, bcf_hdr_t* hts_hdr)
    {
      destination = block();
      if (bcf_read(hts_file, hts_hdr, destination.hts_rec_) >= 0)
      {
        bcf_unpack(destination.hts_rec_, BCF_UN_ALL); //BCF_UN_STR | BCF_UN_FMT);
        bcf_get_genotypes(hts_hdr, destination.hts_rec_, &(destination.gt_), &(destination.gt_sz_));
        destination.num_samples_ = hts_hdr->n[BCF_DT_SAMPLE];
        if (destination.gt_sz_ % destination.num_samples_ != 0)
        {
          // TODO: mixed ploidy at site error.
        }
        else
        {
          std::uint16_t i = 1;
          do
          {
            destination.markers_.emplace_back(destination.hts_rec_, destination.gt_, destination.gt_sz_, i);
            ++i;
          } while (i < destination.hts_rec_->n_allele);
          return true;
        }
      }
      return false;
    }

    reader::reader(const std::string& file_path) :
      hts_file_(bcf_open(file_path.c_str(), "r")), hts_hdr_(nullptr)
    {
      if (hts_file_)
      {
        hts_hdr_ = bcf_hdr_read(hts_file_);
      }
    }

    reader::~reader()
    {
      try
      {
        if (hts_hdr_)
          bcf_hdr_destroy(hts_hdr_);
        if (hts_file_)
          bcf_close(hts_file_);
      }
      catch (...)
      {
        assert(!"bcf_hdr_destroy or bcf_close threw exception!");
      }
    }

    bool reader::read_next_block(block& destination)
    {
      return block::read_block(destination, hts_file_, hts_hdr_);
    }

    char** reader::samples_begin() const
    {
      return hts_hdr_->samples;
    }

    char** reader::samples_end() const
    {
      return hts_hdr_->samples + bcf_hdr_nsamples(hts_hdr_);
    }

    std::string reader::get_chromosome(const reader& rdr, const marker& mkr)
    {
      std::string ret(bcf_hdr_id2name(rdr.hts_hdr_, mkr.chrom_id()));
      return ret;
    }
  }
}
