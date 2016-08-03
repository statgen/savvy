
#include "vcf_reader.hpp"

#include <assert.h>
#include <vcf.h>

namespace vc
{
  namespace vcf
  {
    marker::marker(int* gt, int num_gt, std::uint16_t allele_index) :
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

    const allele_status& marker::operator[](std::size_t i) const
    {
      if (gt_[i] == bcf_gt_missing)
        return const_is_missing;
      else
        return (bcf_gt_allele(gt_[i]) == allele_index_ ? const_has_alt : const_has_ref);
    }

    block::block() :
      hts_rec_(bcf_init1()),
      num_gt_(0),
      gt_(nullptr),
      gt_sz_(0)
    {

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
      if (bcf_read(hts_file, hts_hdr, destination.hts_rec_) >= 0)
      {
        destination.num_gt_ = bcf_get_genotypes(hts_hdr, destination.hts_rec_, &(destination.gt_), &(destination.gt_sz_));
        destination.num_samples_ = hts_hdr->n[BCF_DT_SAMPLE];
        if (destination.num_gt_ % destination.num_samples_ != 0)
        {
          // TODO: mixed ploidy at site error.
        }
        else
        {
          for (std::uint16_t i = 1; i < destination.hts_rec_->n_allele; ++i) // TODO: figure out if n_alleles includes ref.
            destination.markers_.emplace_back(destination.gt_, destination.num_gt_, i);
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
        hts_hdr_ = vcf_hdr_read(hts_file_);
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
  }
}
