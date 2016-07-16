#ifndef VC_M3VCF_READER_HPP
#define VC_M3VCF_READER_HPP

#include <string>
#include <cstdint>
#include <vector>

namespace vc
{


  class m3vcf_block
  {
  private:
    static const char has_alt_code = '1';
    static const char has_ref_code = '0';
    static const char is_missing_code = '.';
  public:
    class marker_iterator
    {
    public:
      class variant_iterator
      {
      public:
        variant_iterator(marker_iterator& parent, std::uint64_t offset);

        bool has_alt_at(std::uint8_t allele_off) const
        {
          return parent_.has_alt_at(offset_, allele_off);
        }

        bool has_ref_at(std::uint8_t allele_off) const
        {
          return parent_.has_ref_at(offset_, allele_off);
        }

        bool is_missing_at(std::uint8_t allele_off) const
        {
          return parent_.is_missing_at(offset_, allele_off);
        }
      private:
        marker_iterator& parent_;
        std::uint64_t offset_;
      };

      marker_iterator(m3vcf_block& parent, std::uint32_t offset);

      bool has_alt_at(std::uint64_t sample_off, std::uint8_t allele_off) const
      {
        return parent_.has_alt_at(offset_, sample_off, allele_off);
      }

      bool has_ref_at(std::uint64_t sample_off, std::uint8_t allele_off) const
      {
        return parent_.has_ref_at(offset_, sample_off, allele_off);
      }

      bool is_missing_at(std::uint64_t sample_off, std::uint8_t allele_off) const
      {
        return parent_.is_missing_at(offset_, sample_off, allele_off);
      }
    private:
      m3vcf_block& parent_;
      std::uint32_t offset_;
    };


    bool has_alt_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t allele_off) const
    {
      return (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[sample_off * ploidy_level_ + allele_off]] == m3vcf_block::has_alt_code);
    }

    bool has_ref_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t allele_off) const
    {
      return (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[sample_off * ploidy_level_ + allele_off]] == m3vcf_block::has_ref_code);
    }

    bool is_missing_at(std::uint32_t marker_off, std::uint64_t sample_off, std::uint8_t allele_off) const
    {
      return (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[sample_off * ploidy_level_ + allele_off]] == m3vcf_block::is_missing_code);
    }


    std::uint64_t sample_size() const;
    std::uint32_t marker_size() const;
  private:
    struct site_info
    {
      std::string chomosome;
      std::uint64_t position;
      std::string ref;
      std::string alt;
    };
    std::vector<site_info> sites_;

    //---- GT Data ----//
    std::vector<std::uint64_t> haplotype_weights_;
    std::vector<std::uint32_t> sample_mappings_;
    std::vector<char> unique_haplotype_matrix_;
    std::uint64_t sample_cnt_;
    std::uint32_t variant_cnt_;
    std::uint32_t unique_haplotype_cnt_;
    std::uint8_t ploidy_level_;
    //---- GT Data ----//
  };

  class m3vcf_reader
  {
  public:
    m3vcf_reader(const std::string& file_path);
    bool read_next_block(m3vcf_block& destination);
  private:
    const std::string file_path_;
  };
}

#endif //VC_M3VCF_READER_HPP
