#include "m3vcf_reader.hpp"

#include <cmath>
#include <assert.h>
#include <cstring>

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
      return parent_.sample_haplotype_at(offset_, i);
    }

    const allele_status& marker::at(std::size_t i) const
    {
      if (i >=  this->haplotype_count())
        throw std::out_of_range("index out of range");
      return (*this)[i];
    }

    std::uint64_t marker::haplotype_count() const
    {
      return parent_.haplotype_count();
    }

    std::uint64_t marker::pos() const
    {
      return position_;
    }

    std::string marker::ref() const
    {
      return ref_;
    }

    std::string marker::alt() const
    {
      return alt_;
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

    const allele_status& block::sample_haplotype_at(std::uint32_t marker_off, std::uint64_t haplotype_off)
    {
      switch (unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[haplotype_off]])
      {
      case '0': return const_has_ref;
      case '1': return const_has_alt;
      case '.': return const_is_missing;
      };
      return const_is_missing;
    }

    const allele_status& block::unique_haplotype_at(std::uint32_t marker_offset, std::uint64_t unique_haplotype_offset)
    {
      switch (unique_haplotype_matrix_[(marker_offset * unique_haplotype_cnt_) + unique_haplotype_offset])
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

    const marker& block::at(std::size_t i) const
    {
      return markers_.at(i);
    }

    bool block::read(block& destination, std::istream& source, std::uint64_t sample_count, std::uint8_t ploidy)
    {
      destination = block();
      std::uint32_t row_length = 0;
      std::uint32_t column_length = 0;
      source.read((char*)&row_length, sizeof(row_length));
      source.read((char*)&column_length, sizeof(column_length));
      row_length = be32toh(row_length);
      column_length = be32toh(column_length);

      destination.unique_haplotype_cnt_ = column_length;


      std::uint8_t byte_width_needed = 4;
      if (destination.unique_haplotype_cnt_ <= 0x100)
        byte_width_needed = 1;
      else if (destination.unique_haplotype_cnt_ <= 0x10000)
        byte_width_needed = 2;

      destination.sample_size_ = sample_count;
      destination.ploidy_level_ = ploidy;

      std::vector<char> buff(destination.haplotype_count() * byte_width_needed);
      destination.sample_mappings_.resize(buff.size(), 0xFFFFFFFF);
      source.read(buff.data(), buff.size());
      switch (byte_width_needed)
      {
        case 1:
        {
          for (std::size_t i = 0; i < destination.haplotype_count(); ++i)
          {
            std::uint8_t byte;
            std::memcpy(&byte, &buff[i], 1);
            destination.sample_mappings_[i] = byte;
          }
          break;
        }
        case 2:
        {
          for (std::size_t i = 0; i < destination.haplotype_count(); ++i)
          {
            std::uint16_t nbo16;
            std::memcpy(&nbo16, &buff[i * 2], 2);
            destination.sample_mappings_[i] = ntohs(nbo16);
          }
          break;
        }
        default:
        {
          for (std::size_t i = 0; i < destination.haplotype_count(); ++i)
          {
            std::uint32_t nbo32;
            std::memcpy(&nbo32, &buff[i * 4], 4);
            destination.sample_mappings_[i] = be32toh(nbo32);
          }
        }
      }

      destination.markers_.reserve(row_length);
      destination.unique_haplotype_matrix_.resize(std::size_t(row_length) * std::size_t(column_length));

      for (std::size_t i = 0; i < row_length; ++i)
      {
        std::uint64_t pos;
        source.read((char*)(&pos), 8);

        std::uint8_t sz = 0;
        std::string marker_id;
        source.read((char*)&sz, 1);
        if (sz)
        {
          marker_id.resize(sz);
          source.read(&marker_id[0], sz);
        }


        sz = 0;
        std::string ref;
        source.read((char*)&sz, 1);
        if (sz)
        {
          ref.resize(sz);
          source.read(&ref[0], sz);
        }

        sz = 0;
        std::string alt;
        source.read((char*)&sz, 1);
        if (sz)
        {
          alt.resize(sz);
          source.read(&alt[0], sz);
        }

        destination.markers_.emplace_back(destination, i, "[CHROM]", be64toh(pos), ref, alt);

        source.read(&(destination.unique_haplotype_matrix_[i * column_length]), column_length);
      }

      return source.good();
    }

    bool block::write(std::ostream& destination, const block& source)
    {

      std::uint32_t row_length;
      std::uint32_t column_length;
      assert(source.markers_.size() <= 0xFFFFFFFF);
      row_length = htobe32(source.markers_.size());
      column_length = htobe32(source.unique_haplotype_cnt_);
      destination.write((char*)&row_length, sizeof(row_length));
      destination.write((char*)&column_length, sizeof(column_length));

      //long byte_width_needed = static_cast<long>(std::ceil(std::log2(source.unique_haplotype_cnt_ + 1) / 8));
      std::uint8_t byte_width_needed = 4;
      if (source.unique_haplotype_cnt_ <= 0x100)
        byte_width_needed = 1;
      else if (source.unique_haplotype_cnt_ <= 0x10000)
        byte_width_needed = 2;

      std::vector<char> buff(source.haplotype_count() * byte_width_needed);
      switch (byte_width_needed)
      {
        case 1:
        {
          for (std::size_t i = 0; i < source.haplotype_count(); ++i)
          {
            std::uint8_t byte = std::uint8_t(source.sample_mappings_[i]);
            std::memcpy(&buff[i], &byte, 1);
          }
          break;
        }
        case 2:
        {
          for (std::size_t i = 0; i < source.haplotype_count(); ++i)
          {
            std::uint16_t nbo16 = htons(source.sample_mappings_[i]);
            std::memcpy(&buff[i * 2], &nbo16, 2);
          }
          break;
        }
        default:
        {
          for (std::size_t i = 0; i < source.haplotype_count(); ++i)
          {
            std::uint32_t nbo32 = htobe32(source.sample_mappings_[i]);
            std::memcpy(&buff[i * 4], &nbo32, 4);
          }
        }
      }

      destination.write(buff.data(), buff.size());

      std::size_t i = 0;
      for (auto it = source.markers_.begin(); it != source.markers_.end(); ++it,++i)
      {
        std::uint64_t pos_nbo = htobe64(it->pos());
        destination.write((char*)(&pos_nbo), 8);
        destination.put(0); // rsid
        std::uint8_t sz = std::uint8_t(0xFF & it->ref().size());
        destination.put(sz);
        destination.write(it->ref().data(), sz);

        sz = std::uint8_t(0xFF & it->alt().size());
        destination.put(sz);
        destination.write(it->alt().data(), sz);

        destination.write(&(source.unique_haplotype_matrix_[i * source.unique_haplotype_cnt_]), source.unique_haplotype_cnt_);
      }


      return destination.good();
    }
    //================================================================//

    //================================================================//
    reader::reader(std::istream& input_stream)
      :
      input_stream_(input_stream)
    {
      std::string version_string(9, '\0');
      input_stream_.read(&version_string[0], version_string.size());

      std::uint8_t str_sz = 0;
      std::string chromosome;
      input_stream_.read((char*)&str_sz, 1);
      if (str_sz)
      {
        chromosome.resize(str_sz);
        input_stream_.read(&chromosome[0], str_sz);
      }


      input_stream_.read((char*)&ploidy_level_, 1);

      std::uint64_t sample_size = 0;
      input_stream_.read((char*)(&sample_size), 8);
      sample_size = be64toh(sample_size);
      sample_ids_.reserve(sample_size);
      for (std::size_t i = 0; i < sample_size; ++i)
      {
        std::uint8_t sz = 0;
        input_stream_.read((char*)(&sz), 1);
        sample_ids_.emplace_back(sz, '\0');
        if (sz)
        {
          input_stream_.read(&sample_ids_.back()[0], sz);
        }
      }
    }

    bool reader::read_next_block(block& destination)
    {
      return block::read(destination, input_stream_, sample_ids_.size(), ploidy_level_);
    }
    //================================================================//
  }
}