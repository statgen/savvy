#include "cvcf_reader.hpp"

namespace vc
{
  namespace cvcf
  {
    //================================================================//
    void marker::for_each_allele(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
    {
      std::uint64_t haplotype_offset = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        haplotype_offset += it->offset;
        if (it->is_allele)
        {
          fn(haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
        }
      }
    }

    void marker::for_each_missing(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
    {
      std::uint64_t haplotype_offset = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        haplotype_offset += it->offset;
        if (!(it->is_allele))
        {
          fn(haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
        }
      }
    }

    void marker::for_each_non_ref(const std::function<void(allele_status status, std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
    {
      std::uint64_t haplotype_offset = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        haplotype_offset += it->offset;
        fn(it->is_allele ? allele_status::has_alt : allele_status::is_missing, haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
      }
    }

    double marker::calculate_allele_frequency() const
    {
      std::uint64_t allele_cnt = 0;
      std::uint64_t total_haplotypes = sample_count_ * ploidy_level_;
      std::uint64_t haplotype_offset = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        if (it->is_allele)
          ++allele_cnt;
        else // missing
          --total_haplotypes;
      }

      return static_cast<double>(allele_cnt) / static_cast<double>(total_haplotypes);
    }
    //================================================================//
  }
}

