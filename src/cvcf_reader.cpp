#include "cvcf_reader.hpp"
#include "varint.hpp"


namespace vc
{
  namespace cvcf
  {
    //================================================================//
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
      std::istreambuf_iterator<char> end_it;
      varint_decode(std::istreambuf_iterator<char>(is), end_it, destination.position_);

      std::uint64_t sz;
      varint_decode(std::istreambuf_iterator<char>(is), end_it, sz);
      destination.id_.resize(sz);
      if (destination.id_.size())
        is.read(&destination.id_[0], destination.alt_.size());

      varint_decode(std::istreambuf_iterator<char>(is), end_it, sz);
      destination.ref_.resize(sz);
      if (destination.ref_.size())
        is.read(&destination.ref_[0], destination.ref_.size());

      varint_decode(std::istreambuf_iterator<char>(is), end_it, sz);
      destination.alt_.resize(sz);
      if (destination.alt_.size())
        is.read(&destination.alt_[0], destination.alt_.size());

      std::uint64_t total_offset = 0;
      destination.non_zero_haplotypes_.reserve(sz);
      for (auto it = destination.non_zero_haplotypes_.begin(); it != destination.non_zero_haplotypes_.end(); ++it,++total_offset)
      {
        std::uint8_t allele;
        std::uint64_t offset;
        one_bit_prefixed_varint::decode(std::istreambuf_iterator<char>(is), end_it, allele, offset);
        total_offset += offset;
        it->offset = total_offset;
        it->status = (allele ? allele_status::has_alt : allele_status::is_missing);
      }

      return is.good();
    }
    //================================================================//
  }
}

