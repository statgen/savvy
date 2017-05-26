#include "savvy/sav_reader.hpp"
#include "savvy/variant_iterator.hpp"

#include <assert.h>
#include <algorithm>


namespace savvy
{
  namespace sav
  {
    //================================================================//
    reader_base::reader_base(const std::string& file_path) :
      input_stream_(file_path),
      file_path_(file_path)
    {
      std::string version_string(7, '\0');
      input_stream_.read(&version_string[0], version_string.size());


      std::istreambuf_iterator<char> in_it(input_stream_);
      std::istreambuf_iterator<char> end;

      std::uint64_t flags{};
      if (varint_decode(in_it, end, flags) != end)
      {
        ++in_it;
        value_bit_width_ = std::uint8_t(flags & 0x07);

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
      }

      input_stream_.peek();
    }

    reader_base::reader_base(reader_base&& source) :
      sample_ids_(std::move(source.sample_ids_)),
      chromosome_(std::move(source.chromosome_)),
      //sbuf_(std::move(source.sbuf_)),
      //input_stream_(&sbuf_),
      input_stream_(std::move(source.input_stream_)),
      file_path_(std::move(source.file_path_)),
      ploidy_level_(source.ploidy_level_),
      metadata_fields_(std::move(source.metadata_fields_)),
      value_bit_width_(source.value_bit_width_)
    {
    }

    reader_base& reader_base::operator=(reader_base&& source)
    {
      if (&source != this)
      {
        sample_ids_ = std::move(source.sample_ids_);
        chromosome_ = std::move(source.chromosome_);
        //sbuf_ = std::move(source.sbuf_);
        //input_stream_.rdbuf(&sbuf_);
        input_stream_ = std::move(source.input_stream_);
        file_path_ = std::move(source.file_path_);
        ploidy_level_ = source.ploidy_level_;
        metadata_fields_ = std::move(source.metadata_fields_);
        value_bit_width_ = source.value_bit_width_;
      }
      return *this;
    }

//    reader& reader_base::operator>>(marker& destination)
//    {
//      marker::read(destination, sample_ids_.size() * ploidy_level_, input_stream_);
//      return *this;
//    }
    //================================================================//

    //================================================================//
    const std::array<std::string, 0> writer::empty_string_array = {};

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


      variant_iterator<std::vector<float>> it(r);
      variant_iterator<std::vector<float>> end;
      while (it != end && start_pos >= 0)
      {
        std::int64_t end_pos = r.tellg();
        if (start_pos >= 0 && end_pos >= 0)
        {
          s1r::index_base::entry e(it->locus(), it->locus() + std::max(it->ref().size(), it->alt().size()) - 1, static_cast<std::uint64_t>(start_pos), static_cast<std::uint64_t>(end_pos) - static_cast<std::uint64_t>(start_pos));
          max_region_value = std::max(max_region_value, e.region_end());
          max_file_position = std::max(max_file_position, e.value().first);
          max_file_position = std::max(max_file_position, e.value().second);
          entries.emplace_back(std::move(e));
        }
        start_pos = end_pos;
        ++i;
        ++it;
      }

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

