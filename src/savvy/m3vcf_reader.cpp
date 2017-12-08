/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/m3vcf_reader.hpp"

#include <cmath>
#include <assert.h>
#include <cstring>

namespace savvy
{
  namespace m3vcf
  {
    namespace detail
    {
      std::uint64_t ceil_divide(std::uint64_t dividend, std::uint64_t divisor)
      {
        return (std::uint64_t)(1) + ((dividend - (std::uint64_t)(1)) / divisor);
      }

      void deserialize_string(std::string& output, std::istream& input)
      {
        std::size_t i = 0;
        char cur = 0;
        while (input.get(cur) && cur == (char)255)
          ++i;
        output.resize(i * 255 + (unsigned char)cur);
        if (output.size())
          input.read(&output[0], output.size());
      }

      void serialize_string(std::ostream& output, const std::string& input)
      {
        std::size_t string_size = input.size();
        while (string_size >= 255)
        {
          output.put((char)255);
          string_size -= 255;
        }
        output.put((char)string_size);
        output.write(input.data(), input.size());
      }
    }

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
      switch (unique_haplotype_matrix_[sample_mappings_[haplotype_off]][marker_off]) //unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + sample_mappings_[haplotype_off]])
      {
      case '0': return const_has_ref;
      case '1': return const_has_alt;
      case '.': return const_is_missing;
      };
      return const_is_missing;
    }

    const allele_status& block::unique_haplotype_at(std::uint32_t marker_offset, std::uint64_t unique_haplotype_offset)
    {
      switch (unique_haplotype_matrix_[unique_haplotype_offset][marker_offset]) //unique_haplotype_matrix_[(marker_offset * unique_haplotype_cnt_) + unique_haplotype_offset])
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
        char hap =  unique_haplotype_matrix_[i][marker_off]; //unique_haplotype_matrix_[(marker_off * unique_haplotype_cnt_) + i];
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
      std::uint32_t num_rows = 0;
      std::uint32_t num_columns = 0;
      source.read((char*)&num_rows, sizeof(num_rows));
      source.read((char*)&num_columns, sizeof(num_columns));
      num_rows = be32toh(num_rows);
      num_columns = be32toh(num_columns);

      destination.unique_haplotype_cnt_ = num_columns;


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
            destination.sample_mappings_[i] = be16toh(nbo16);
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

      destination.markers_.reserve(num_rows);
      //destination.unique_haplotype_matrix_.resize(std::size_t(num_rows) * std::size_t(num_columns));
      destination.unique_haplotype_matrix_.resize(std::size_t(num_columns), std::vector<char>(std::size_t(num_rows), '\0'));

      for (std::size_t i = 0; i < num_rows; ++i)
      {
        std::uint64_t pos;
        source.read((char*)(&pos), 8);


        std::string marker_id;
        detail::deserialize_string(marker_id, source);

        std::string ref;
        detail::deserialize_string(ref, source);

        std::string alt;
        detail::deserialize_string(alt, source);

        destination.markers_.emplace_back(destination, i, "[CHROM]", be64toh(pos), ref, alt);

        //source.read(&(destination.unique_haplotype_matrix_[i * num_columns]), num_columns);
        std::vector<char> tmp(num_columns);
        source.read(tmp.data(), num_columns);
        for (std::size_t j = 0; j < num_columns; ++j)
          destination.unique_haplotype_matrix_[j][i] = tmp[j];

      }

      return source.good();
    }

    bool block::write(std::ostream& destination, const block& source)
    {

      std::uint32_t num_rows;
      std::uint32_t num_columns;
      assert(source.markers_.size() <= 0xFFFFFFFF);
      num_rows = htobe32(source.markers_.size());
      num_columns = htobe32(source.unique_haplotype_cnt_);
      destination.write((char*)&num_rows, sizeof(num_rows));
      destination.write((char*)&num_columns, sizeof(num_columns));

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
            std::uint16_t nbo16 = htobe16(source.sample_mappings_[i]);
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

        detail::serialize_string(destination, ""); // rsid

        detail::serialize_string(destination, it->ref());

        detail::serialize_string(destination, it->alt());

        //destination.write(&(source.unique_haplotype_matrix_[i * source.unique_haplotype_cnt_]), source.unique_haplotype_cnt_);
        std::vector<char> tmp(source.unique_haplotype_matrix_.size());
        for (std::size_t j = 0; j < source.unique_haplotype_matrix_.size(); ++j)
          tmp[j] = source.unique_haplotype_matrix_[j][i];
        destination.write(tmp.data(), tmp.size());
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


      std::string chromosome;
      detail::deserialize_string(chromosome, input_stream_);

      input_stream_.read((char*)(&ploidy_level_), 4);
      ploidy_level_ = be32toh(ploidy_level_);

      std::uint32_t sample_size = 0;
      input_stream_.read((char*)(&sample_size), 4);
      sample_size = be32toh(sample_size);
      sample_ids_.reserve(sample_size);
      for (std::size_t i = 0; i < sample_size; ++i)
      {
        sample_ids_.emplace_back();
        detail::deserialize_string(sample_ids_.back(), input_stream_);
      }
    }

    reader& reader::operator>>(block& destination)
    {
      block::read(destination, input_stream_, sample_ids_.size(), ploidy_level_);
      return  *this;
    }
    //================================================================//
  }
}