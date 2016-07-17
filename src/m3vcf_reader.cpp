#include "m3vcf_reader.hpp"

#include <cmath>

std::uint64_t ceil_divide(std::uint64_t dividend, std::uint64_t divisor)
{
  return (std::uint64_t)(1) + ((dividend - (std::uint64_t)(1)) / divisor);
}

namespace vc
{
  namespace m3vcf
  {
    //================================================================//
    marker::marker(block& parent, std::uint32_t offset, const std::string& chromosome, std::uint64_t position, const std::string& ref, const std::string& alt)
      :
      parent_(parent),
      offset_(offset),
      chromosome_(chromosome),
      position_(position),
      ref_(ref),
      alt_(alt)
    {
    }

    bool marker::has_alt_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return parent_.has_alt_at(offset_, sample_off, ploidy_off);
    }

    bool marker::has_ref_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return parent_.has_ref_at(offset_, sample_off, ploidy_off);
    }

    bool marker::is_missing_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return parent_.is_missing_at(offset_, sample_off, ploidy_off);
    }

    allele_status marker::operator()(std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return parent_(offset_, sample_off, ploidy_off);
    }
    //================================================================//

    //================================================================//
    bool block::has_alt_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return (operator()(marker_off, sample_off, ploidy_off) == allele_status::is_missing);
    }

    bool block::has_ref_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return (operator()(marker_off, sample_off, ploidy_off) == allele_status::is_missing);
    }

    bool block::is_missing_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t ploidy_off) const
    {
      return (operator()(marker_off, sample_off, ploidy_off) == allele_status::is_missing);
    }

    allele_status block::operator()(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t allele_off) const
    {
      switch (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[sample_off * ploidy_level_ + allele_off]])
      {
      case '0': return allele_status::has_ref;
      case '1': return allele_status::has_alt;
      case '.': return allele_status::is_missing;
      };
      return allele_status::is_missing;
    }

    block::const_iterator block::begin()
    {
      return const_iterator(this->markers_.data());
    }

    block::const_iterator block::end()
    {
      return const_iterator(this->markers_.data() + this->markers_.size());
    }
    //================================================================//

    //================================================================//
    reader::reader(const std::string& file_path)
      :
      file_path_(file_path)
    {

    }
    //================================================================//
  }
}