#include "cvcf_reader.hpp"
#include "varint.hpp"

#include <iostream>

namespace vc
{
  namespace cvcf
  {
    //================================================================//
//    void marker::for_each_allele(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
//    {
//      std::uint64_t haplotype_offset = 0;
//      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
//      {
//        haplotype_offset += it->offset;
//        if (it->is_allele)
//        {
//          fn(haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
//        }
//      }
//    }
//
//    void marker::for_each_missing(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
//    {
//      std::uint64_t haplotype_offset = 0;
//      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
//      {
//        haplotype_offset += it->offset;
//        if (!(it->is_allele))
//        {
//          fn(haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
//        }
//      }
//    }
//
//    void marker::for_each_non_ref(const std::function<void(allele_status status, std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn)
//    {
//      std::uint64_t haplotype_offset = 0;
//      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
//      {
//        haplotype_offset += it->offset;
//        fn(it->is_allele ? allele_status::has_alt : allele_status::is_missing, haplotype_offset / ploidy_level_, (std::uint8_t)(haplotype_offset % ploidy_level_));
//      }
//    }

    const marker::const_iterator::value_type marker::const_iterator::const_is_missing = allele_status::is_missing;
    const marker::const_iterator::value_type marker::const_iterator::const_has_ref = allele_status::has_ref;
    const marker::const_iterator::value_type marker::const_iterator::const_has_alt = allele_status::has_alt;

    marker::non_ref_iterator marker::non_ref_begin() const
    {
      return non_zero_haplotypes_.begin();
    }

    marker::non_ref_iterator marker::non_ref_end() const
    {
      return non_zero_haplotypes_.end();
    }

    marker::const_iterator marker::begin() const
    {
      return const_iterator(0, non_zero_haplotypes_.begin(), non_zero_haplotypes_.end());
    }

    marker::const_iterator marker::end() const
    {
      return const_iterator(sample_count_ * ploidy_level_, non_zero_haplotypes_.begin(), non_zero_haplotypes_.end());
    }

    double marker::calculate_allele_frequency() const
    {
      std::uint64_t allele_cnt = 0;
      std::uint64_t total_haplotypes = sample_count_ * ploidy_level_;
      std::uint64_t haplotype_offset = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        if (it->status == allele_status::is_missing)
          --total_haplotypes;
        else // has alt
          ++allele_cnt;
      }

      return static_cast<double>(allele_cnt) / static_cast<double>(total_haplotypes);
    }

    bool marker::read(marker& destination, std::istream& is)
    {
      std::istreambuf_iterator<char> input_it(is);
      destination.position_ = varint_decode(input_it);

      destination.id_.resize(varint_decode(input_it));
      if (destination.id_.size())
        is.read(&destination.id_[0], destination.alt_.size());

      destination.ref_.resize(varint_decode(input_it));
      if (destination.ref_.size())
        is.read(&destination.ref_[0], destination.ref_.size());

      destination.alt_.resize(varint_decode(input_it));
      if (destination.alt_.size())
        is.read(&destination.alt_[0], destination.alt_.size());

      std::uint64_t offset = 0;
      destination.non_zero_haplotypes_.reserve(varint_decode(input_it));
      for (auto it = destination.non_zero_haplotypes_.begin(); it != destination.non_zero_haplotypes_.end(); ++it,++offset)
      {
        std::uint8_t allele = 0;
        offset += one_bit_prefixed_varint::decode(allele, input_it);
        it->offset = offset;
        it->status = (allele ? allele_status::has_alt : allele_status::is_missing);
      }

      return true;
    }
    //================================================================//
  }
}

