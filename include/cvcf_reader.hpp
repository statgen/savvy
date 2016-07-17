#ifndef LIBVC_CVCF_READER_HPP
#define LIBVC_CVCF_READER_HPP

#include "allele_status.hpp"

#include <cstdint>
#include <string>
#include <vector>
#include <functional>

namespace vc
{
  namespace cvcf
  {
    class marker
    {
    public:
      bool has_alt_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
      bool has_ref_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
      bool is_missing_at(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
      allele_status operator()(std::uint64_t sample_off, std::uint8_t ploidy_off) const;
      void for_each_allele(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
      void for_each_missing(const std::function<void(std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
      void for_each_non_ref(const std::function<void(allele_status status, std::uint64_t sample_off, std::uint8_t ploidy_off)>& fn);
      double calculate_allele_frequency() const;
    private:
      struct hap_pair
      {
        std::uint64_t offset;
        bool is_allele;
      };
      std::vector<hap_pair> non_zero_haplotypes_;
      std::uint8_t ploidy_level_;
      std::uint64_t sample_count_;
    };

    class reader
    {
    public:
      bool read_next_marker(marker& destination);
    private:
      std::uint8_t ploidy_level_;
      std::uint64_t sample_count_;
    };
  }
}

#endif //LIBVC_CVCF_READER_HPP
