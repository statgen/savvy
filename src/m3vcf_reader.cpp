#include "m3vcf_reader.hpp"

#include <unordered_set>
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

    const allele_status& marker::operator[](std::size_t i) const
    {
      return parent_.haplotype_at(offset_, i);
    }

    std::uint64_t marker::haplotype_count() const
    {
      return parent_.haplotype_count();
    }

    double marker::calculate_allele_frequency() const
    {
      return parent_.calculate_allele_frequency(offset_);
    }
    //================================================================//

    //================================================================//
    const allele_status block::const_has_ref = allele_status::has_ref;
    const allele_status block::const_has_alt = allele_status::has_alt;
    const allele_status block::const_is_missing = allele_status::is_missing;

    const allele_status& block::haplotype_at(std::uint32_t marker_off, std::uint64_t haplotype_off)
    {
      switch (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[haplotype_off]])
      {
      case '0': return const_has_ref;
      case '1': return const_has_alt;
      case '.': return const_is_missing;
      };
      return const_is_missing;
    }

    double block::calculate_allele_frequency(std::uint32_t marker_off) const
    {
      std::uint64_t allele_cnt = 0;
      std::uint64_t total_haplotypes = sample_size_ * ploidy_level_;

      for (std::uint32_t i = 0; i < unique_haplotype_cnt_; ++i)
      {
        char hap = unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + i];
        if (hap == '1')
          allele_cnt += haplotype_weights_[i];
        else if (hap != '0') // missing
          total_haplotypes -= haplotype_weights_[i];
      }

      return static_cast<double>(allele_cnt) / static_cast<double>(total_haplotypes);
    }

    block::const_iterator block::begin()
    {
      return const_iterator(this->markers_.data());
    }

    block::const_iterator block::end()
    {
      return const_iterator(this->markers_.data() + this->markers_.size());
    }

    const marker& block::operator[](std::size_t i) const
    {
      return markers_[i];
    }

    bool block::add_marker(std::uint64_t position, const std::string& ref, const std::string& alt, const char* hap_array, std::size_t hap_array_sz)
    {
      bool ret = false;

      if (hap_array_sz == sample_size_ * ploidy_level_)
      {
        std::int64_t current_savings = static_cast<std::int64_t>(hap_array_sz * markers_.size()) - static_cast<std::int64_t>(unique_haplotype_cnt_ * markers_.size() + hap_array_sz);

        std::unordered_set<std::string> new_unique_haps;
        std::string column_string(markers_.size() + 1, '\0');

        for (std::size_t i = 0; i < hap_array_sz; ++i)
        {
          std::size_t j = 0;
          for (; j < markers_.size(); ++j)
          {
            column_string[j] = unique_haplotype_matrix_[(j * unique_haplotype_cnt_) + sample_mappings_[i]];
          }
          column_string[j] = hap_array[i];
          new_unique_haps.insert(column_string);
        }

        std::int64_t new_savings = static_cast<std::int64_t>(hap_array_sz * (markers_.size() + 1)) - static_cast<std::int64_t>(new_unique_haps.size() * (markers_.size() + 1) + hap_array_sz);
        if (new_savings >= current_savings)
        {
          markers_.push_back(marker(*this, markers_.size(), "", position, ref, alt));
          std::vector<char> new_unqiue_haplotype_matrix(new_unique_haps.size() * markers_.size());
          std::size_t i = 0;
          for (auto it = new_unique_haps.begin(); it != new_unique_haps.end(); ++it,++i)
          {
            for (std::size_t j = 0; j < it->size(); ++j)
            {
              new_unqiue_haplotype_matrix[(j * new_unique_haps.size()) + i] = (*it)[j];
            }
          }

          unique_haplotype_matrix_ = std::move(new_unqiue_haplotype_matrix);
          ret = true;
        }
      }

      return ret;
    }
    //================================================================//

    //================================================================//
    reader::reader(std::istream& input_stream)
      :
      input_stream_(input_stream)
    {

    }

    bool reader::read_next_block(block& destination)
    {
      return block::read_block(destination, input_stream_);
    }
    //================================================================//
  }
}