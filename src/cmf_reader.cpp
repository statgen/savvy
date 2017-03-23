#include "cmf_reader.hpp"

#include <assert.h>
#include <algorithm>



namespace vc
{
  namespace cmf
  {
    //================================================================//
    const marker::const_iterator::value_type marker::const_iterator::const_is_missing = allele_status::is_missing;
    const marker::const_iterator::value_type marker::const_iterator::const_has_ref = allele_status::has_ref;
    const marker::const_iterator::value_type marker::const_iterator::const_has_alt = allele_status::has_alt;

    const allele_status const_is_missing = allele_status::is_missing;
    const allele_status const_has_ref = allele_status::has_ref;
    const allele_status const_has_alt = allele_status::has_alt;

    const allele_status& marker::operator[](std::uint64_t i) const
    {
      auto end = non_zero_haplotypes_.end();
      sparse_vector_allele val;
      val.offset = i;
      auto it = std::lower_bound(non_zero_haplotypes_.begin(), end, val, [](const sparse_vector_allele& a, const sparse_vector_allele& b) { return a.offset < b.offset; });
      if (it == end || it->offset != i)
        return const_has_ref;
      else
        return it->status;
    }

    const allele_status& marker::at(std::uint64_t i) const
    {
      if (i >= non_zero_haplotypes_.size())
        throw std::out_of_range("index out of range");
      return (*this)[i];
    }

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
      return const_iterator(haplotype_count_, non_zero_haplotypes_.begin(), non_zero_haplotypes_.end());
    }

    double marker::calculate_allele_frequency() const
    {
      std::uint64_t allele_cnt = 0;
      std::uint64_t total_haplotypes = haplotype_count_;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        if (it->status == allele_status::is_missing)
          --total_haplotypes;
        else // has alt
          ++allele_cnt;
      }

      return static_cast<double>(allele_cnt) / static_cast<double>(total_haplotypes);
    }

    void marker::read(marker& destination, std::uint64_t haplotype_count, std::istream& is)
    {
      destination.haplotype_count_ = haplotype_count;

      std::istreambuf_iterator<char> in_it(is);
      std::istreambuf_iterator<char> end_it;

      if (varint_decode(in_it, end_it, destination.position_) != end_it)
      {
        ++in_it;
        std::uint64_t sz;
        if (varint_decode(in_it, end_it, sz) != end_it)
        {
          ++in_it;
          destination.ref_.resize(sz);
          if (sz)
            is.read(&destination.ref_[0], destination.ref_.size());


          if (varint_decode(in_it, end_it, sz) != end_it)
          {
            ++in_it;
            destination.alt_.resize(sz);
            if (sz)
              is.read(&destination.alt_[0], destination.alt_.size());

            // TODO: Read metadata values.

            varint_decode(in_it, end_it, sz);
            destination.non_zero_haplotypes_.resize(sz);

            std::uint64_t total_offset = 0;
            for (auto it = destination.non_zero_haplotypes_.begin(); it != destination.non_zero_haplotypes_.end() && in_it != end_it; ++it,++total_offset)
            {
              std::uint8_t allele;
              std::uint64_t offset;
              one_bit_prefixed_varint::decode(++in_it, end_it, allele, offset);
              total_offset += offset;
              it->offset = total_offset;
              it->status = (allele ? allele_status::has_alt : allele_status::is_missing);
            }

//            std::uint8_t rle = 0;
//            one_bit_prefixed_varint::decode(in_it, end_it, rle, sz);
//            if (rle)
//            {
//              std::vector<std::pair<sparse_vector_allele, std::uint64_t>> rle_non_zero_haplotypes(sz);
//              std::uint64_t total_repeat_bytes = 0;
//              for (auto it = rle_non_zero_haplotypes.begin(); it != rle_non_zero_haplotypes.end() && in_it != end_it; ++it)
//              {
//                std::uint8_t allele_repeat_prefix;
//                two_bit_prefixed_varint::decode(++in_it, end_it, allele_repeat_prefix, it->first.offset);
//                it->first.status = (allele_repeat_prefix & 0x80 ? allele_status::has_alt : allele_status::is_missing);
//
//                if (allele_repeat_prefix & 0x40 && in_it != end_it)
//                {
//                  // Read repeat byte.
//                  varint_decode(++in_it, end_it, it->second);
//                  total_repeat_bytes += it->second;
//                }
//                else
//                {
//                  it->second = 0;
//                }
//              }
//
//              std::uint64_t total_offset = 0;
//              destination.non_zero_haplotypes_.clear();
//              destination.non_zero_haplotypes_.reserve(sz + total_repeat_bytes);
//              for (auto it = rle_non_zero_haplotypes.begin(); it != rle_non_zero_haplotypes.end(); ++it,++total_offset)
//              {
//                total_offset += it->first.offset;
//                destination.non_zero_haplotypes_.emplace_back(it->first.status, total_offset);
//                while (it->second)
//                {
//                  total_offset += (it->first.offset + 1);
//                  destination.non_zero_haplotypes_.emplace_back(it->first.status, total_offset);
//                  --(it->second);
//                }
//              }
//            }
//            else
//            {
//              destination.non_zero_haplotypes_.resize(sz);
//
//              std::uint64_t total_offset = 0;
//              for (auto it = destination.non_zero_haplotypes_.begin(); it != destination.non_zero_haplotypes_.end() && in_it != end_it; ++it,++total_offset)
//              {
//                std::uint8_t allele;
//                std::uint64_t offset;
//                one_bit_prefixed_varint::decode(++in_it, end_it, allele, offset);
//                total_offset += offset;
//                it->offset = total_offset;
//                it->status = (allele ? allele_status::has_alt : allele_status::is_missing);
//              }
//            }
          }
        }
      }

      is.get();
    }

    std::size_t marker::calculate_serialized_gt_size() const
    {
      std::size_t ret = 0;

      std::uint64_t last_pos = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        std::uint64_t offset = it->offset - last_pos;
        last_pos = it->offset + 1;
        ret += one_bit_prefixed_varint::encoded_byte_width(offset);
      }

      return ret;
    }

    std::tuple<std::size_t, std::size_t> marker::calculate_rle_serialized_gt_size_and_count() const
    {
      std::size_t ret_sz = 0;
      std::size_t total_number_of_repeats = 0;

      std::uint64_t last_pos = 0;
      for (auto it = non_zero_haplotypes_.begin(); it != non_zero_haplotypes_.end(); ++it)
      {
        std::uint64_t offset = it->offset - last_pos;
        last_pos = it->offset + 1;

        std::size_t number_of_repeats = 0;
        while (it + 1 != non_zero_haplotypes_.end())
        {
          auto next_it = (it + 1);
          if (offset != next_it->offset - last_pos || it->status != next_it->status)
            break;

          ++number_of_repeats;
          ++it;
          last_pos = it->offset + 1;
        }

        ret_sz += two_bit_prefixed_varint::encoded_byte_width(offset);

        if (number_of_repeats)
        {
          ret_sz += varint_encoded_byte_width(number_of_repeats);
          total_number_of_repeats += number_of_repeats;
        }
      }

      return std::make_tuple(ret_sz, non_zero_haplotypes_.size() - total_number_of_repeats);
    }

    void marker::write(std::ostream& os, const marker& source)
    {
      std::ostreambuf_iterator<char> os_it(os.rdbuf());
      varint_encode(source.position_, os_it);

      varint_encode(source.ref_.size(), os_it);
      if (source.ref_.size())
        std::copy(source.ref_.begin(), source.ref_.end(), os_it);
        //os.write(&source.ref_[0], source.ref_.size());

      varint_encode(source.alt_.size(), os_it);
      if (source.alt_.size())
        std::copy(source.alt_.begin(), source.alt_.end(), os_it);
        //os.write(&source.alt_[0], source.alt_.size());

      varint_encode(source.non_zero_haplotypes_.size(), os_it);
      std::uint64_t last_pos = 0;
      for (auto it = source.non_zero_haplotypes_.begin(); it != source.non_zero_haplotypes_.end(); ++it)
      {
        std::uint64_t offset = it->offset - last_pos;
        last_pos = it->offset + 1;
        std::uint8_t allele = (it->status == allele_status::has_alt ? std::uint8_t(0x80) : std::uint8_t(0x00));
        one_bit_prefixed_varint::encode(allele, offset, os_it);
      }

//      std::uint64_t normal_sz = source.calculate_serialized_gt_size();
//      std::uint64_t rle_size, rle_cnt;
//      std::tie(rle_size, rle_cnt) = source.calculate_rle_serialized_gt_size_and_count();
//      if (rle_size < normal_sz)
//      {
//        std::uint8_t rle = 0x80;
//        one_bit_prefixed_varint::encode(rle, rle_cnt, os_it);
//
//        std::uint64_t last_pos = 0;
//        for (auto it = source.non_zero_haplotypes_.begin(); it != source.non_zero_haplotypes_.end(); ++it)
//        {
//          std::uint64_t offset = it->offset - last_pos;
//          last_pos = it->offset + 1;
//
//          std::size_t number_of_repeats = 0;
//          while (it + 1 != source.non_zero_haplotypes_.end())
//          {
//            auto next_it = (it + 1);
//            if (offset != next_it->offset - last_pos || it->status != next_it->status)
//              break;
//
//            ++number_of_repeats;
//            ++it;
//            last_pos = it->offset + 1;
//          }
//
//          std::uint8_t allele_repeat_prefix = (it->status == allele_status::has_alt ? std::uint8_t(0x80) : std::uint8_t(0x00));
//          if (number_of_repeats)
//            allele_repeat_prefix |= 0x40;
//          two_bit_prefixed_varint::encode(allele_repeat_prefix, offset, os_it);
//
//          if (number_of_repeats)
//            varint_encode(number_of_repeats, os_it);
//        }
//      }
//      else
//      {
//        std::uint8_t rle = 0x00;
//        one_bit_prefixed_varint::encode(rle, source.non_zero_haplotypes_.size(), os_it);
//
//        std::uint64_t last_pos = 0;
//        for (auto it = source.non_zero_haplotypes_.begin(); it != source.non_zero_haplotypes_.end(); ++it)
//        {
//          std::uint64_t offset = it->offset - last_pos;
//          last_pos = it->offset + 1;
//          std::uint8_t allele = (it->status == allele_status::has_alt ? std::uint8_t(0x80) : std::uint8_t(0x00));
//          one_bit_prefixed_varint::encode(allele, offset, os_it);
//        }
//      }
    }
    //================================================================//

    //================================================================//
    reader::reader(const std::string& file_path) :
      sbuf_(file_path),
      input_stream_(&sbuf_),
      file_path_(file_path)
    {
      std::string version_string(7, '\0');
      input_stream_.read(&version_string[0], version_string.size());


      std::istreambuf_iterator<char> in_it(input_stream_);
      std::istreambuf_iterator<char> end;

      std::uint64_t sz;
      if (varint_decode(in_it, end, sz) != end)
      {
        ++in_it;
        if (sz)
        {
          chromosome_.resize(sz);
          input_stream_.read(&chromosome_[0], sz);
        }

        varint_decode(in_it, end, sz);
        assert(sz < 256);
        ploidy_level_ = static_cast<std::uint8_t>(sz);

        if (in_it != end)
        {
          std::uint64_t sample_size;
          if (varint_decode(++in_it, end, sample_size) != end)
          {
            ++in_it;
            sample_ids_.reserve(sample_size);

            std::uint64_t id_sz;
            while (sample_size && varint_decode(in_it, end, id_sz) != end)
            {
              ++in_it;
              sample_ids_.emplace_back();
              if (id_sz)
              {
                sample_ids_.back().resize(id_sz);
                input_stream_.read(&sample_ids_.back()[0], id_sz);
              }
              --sample_size;
            }

            std::uint64_t metadata_fields_cnt;
            if (varint_decode(in_it, end, metadata_fields_cnt) != end)
            {
              ++in_it;
              metadata_fields_.reserve(metadata_fields_cnt);

              std::uint64_t field_sz;
              while (metadata_fields_cnt && varint_decode(in_it, end, field_sz) != end)
              {
                ++in_it;
                metadata_fields_.emplace_back();
                if (field_sz)
                {
                  metadata_fields_.back().resize(field_sz);
                  input_stream_.read(&metadata_fields_.back()[0], field_sz);
                }
                --metadata_fields_cnt;
              }

              if (!metadata_fields_cnt)
                return; //TODO: This is ugly. Consider not depending on on istream error handling.
            }
          }
        }
      }

      input_stream_.peek();
    }

    reader::reader(reader&& source) :
      sample_ids_(std::move(source.sample_ids_)),
      chromosome_(std::move(source.chromosome_)),
      sbuf_(std::move(source.sbuf_)),
      input_stream_(&sbuf_),
      file_path_(std::move(source.file_path_)),
      ploidy_level_(source.ploidy_level_),
      metadata_fields_(std::move(source.metadata_fields_))
    {
    }

    reader& reader::operator=(reader&& source)
    {
      if (&source != this)
      {
        sample_ids_ = std::move(source.sample_ids_);
        chromosome_ = std::move(source.chromosome_);
        sbuf_ = std::move(source.sbuf_);
        input_stream_.rdbuf(&sbuf_);
        file_path_ = std::move(source.file_path_);
        ploidy_level_ = source.ploidy_level_;
        metadata_fields_ = std::move(source.metadata_fields_);
      }
      return *this;
    }

    reader& reader::operator>>(marker& destination)
    {
      marker::read(destination, sample_ids_.size() * ploidy_level_, input_stream_);
      return *this;
    }
    //================================================================//

    //================================================================//
    bool writer::create_index(const std::string& input_file_path, std::string output_file_path)
    {
      bool ret = false;
      std::size_t i = 0;

      if (output_file_path.empty())
        output_file_path = input_file_path + ".s1r";

      std::uint64_t max_region_value = 0;
      std::uint64_t max_file_position = 0;

      std::vector<s1r::index_base::entry> entries;
      reader r(input_file_path);
      std::int64_t start_pos = r.tellg();

      marker buf;
      reader::input_iterator it(r, buf);
      while (it != reader::input_iterator{} && start_pos >= 0)
      {
        std::int64_t end_pos = r.tellg();
        if (start_pos >= 0 && end_pos >= 0)
        {
          s1r::index_base::entry e(it->pos(), it->pos() + std::max(it->ref().size(), it->alt().size()) - 1, static_cast<std::uint64_t>(start_pos), static_cast<std::uint64_t>(end_pos) - static_cast<std::uint64_t>(start_pos));
          max_region_value = std::max(max_region_value, e.region_end());
          max_file_position = std::max(max_file_position, e.value().first);
          max_file_position = std::max(max_file_position, e.value().second);
          entries.emplace_back(std::move(e));
        }
        start_pos = end_pos;
        ++i;
        ++it;
      }
      std::cout << i << std::endl;

      if (start_pos < 0)
      {
        // TODO: handle error.
      }
      else
      {
        ret = s1r::create_file(output_file_path, entries.begin(), entries.end(), s1r::block_size::bs_4096);
      }

      return ret;
    }
    //================================================================//
  }
}

