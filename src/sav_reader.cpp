#include "savvy/sav_reader.hpp"
#include "savvy/variant_iterator.hpp"
#include "savvy/utility.hpp"

#include <assert.h>
#include <algorithm>
#include <map>


namespace savvy
{
  namespace sav
  {
    //================================================================//
    reader_base::reader_base(const std::string& file_path, fmt requested_data_format) :
      input_stream_(file_path),
      file_path_(file_path),
      file_data_format_(fmt::genotype),
      requested_data_formats_(requested_data_format)
    {
      std::string version_string(7, '\0');
      input_stream_.read(&version_string[0], version_string.size());

      std::string uuid(16, '\0');
      input_stream_.read(&uuid[0], uuid.size());


      std::istreambuf_iterator<char> in_it(input_stream_);
      std::istreambuf_iterator<char> end;

      std::uint64_t headers_size;
      if (varint_decode(in_it, end, headers_size) != end)
      {
        ++in_it;
        headers_.reserve(headers_size);

        while (headers_size && in_it != end)
        {
          std::uint64_t key_size;
          if (varint_decode(in_it, end, key_size) != end)
          {
            ++in_it;
            if (key_size)
            {
              std::string key;
              key.resize(key_size);
              input_stream_.read(&key[0], key_size);

              std::uint64_t val_size;
              if (varint_decode(in_it, end, val_size) != end)
              {
                ++in_it;
                if (key_size)
                {
                  std::string val;
                  val.resize(val_size);
                  input_stream_.read(&val[0], val_size);

                  if (key == "INFO")
                  {
                    std::string info_field = parse_header_id(val);
                    metadata_fields_.push_back(std::move(info_field));
                  }
                  else if (key == "FORMAT")
                  {
                    std::string format_field = parse_header_id(val);
                    if (format_field == "GT")
                      file_data_format_ = fmt::genotype;
                    else if (format_field == "GP")
                      file_data_format_ = fmt::genotype_probability;
                  }
                  headers_.emplace_back(std::move(key), std::move(val));
                }
              }

            }
          }
          --headers_size;
        }

        if (!headers_size)
        {
          std::uint64_t sample_size;
          if (varint_decode(in_it, end, sample_size) != end)
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

            if (!sample_size)
              return;
          }
        }
      }

      input_stream_.peek();
    }
#if !defined(__GNUC__) || defined(__clang__) || __GNUC__ > 4
    reader_base::reader_base(reader_base&& source) :
      sample_ids_(std::move(source.sample_ids_)),
      //sbuf_(std::move(source.sbuf_)),
      //input_stream_(&sbuf_),
      input_stream_(std::move(source.input_stream_)),
      file_path_(std::move(source.file_path_)),
      metadata_fields_(std::move(source.metadata_fields_)),
      file_data_format_(source.file_data_format_),
      requested_data_formats_(source.requested_data_formats_)
    {
    }

    reader_base& reader_base::operator=(reader_base&& source)
    {
      if (&source != this)
      {
        sample_ids_ = std::move(source.sample_ids_);
        //sbuf_ = std::move(source.sbuf_);
        //input_stream_.rdbuf(&sbuf_);
        input_stream_ = std::move(source.input_stream_);
        file_path_ = std::move(source.file_path_);
        metadata_fields_ = std::move(source.metadata_fields_);
        file_data_format_ = source.file_data_format_;
        requested_data_formats_ = source.requested_data_formats_;
      }
      return *this;
    }
#endif
//    reader& reader_base::operator>>(marker& destination)
//    {
//      marker::read(destination, sample_ids_.size() * ploidy_level_, input_stream_);
//      return *this;
//    }
    //================================================================//

    //================================================================//
    const std::array<std::string, 0> writer::empty_string_array = {};
    const std::array<std::pair<std::string, std::string>, 0> writer::empty_string_pair_array = {};

    bool writer::create_index(const std::string& input_file_path, std::string output_file_path)
    {
      bool ret = false;

      if (output_file_path.empty())
        output_file_path = input_file_path + ".s1r";

      reader r(input_file_path);
      std::int64_t start_pos = r.tellg();

      std::uint32_t min = std::numeric_limits<std::uint32_t>::max();
      std::uint32_t max = 0;
      std::map<std::string, std::vector<s1r::tree_base::entry>> index_data;
      basic_variant_iterator<reader, dense_allele_vector<float>> it(r);
      basic_variant_iterator<reader, dense_allele_vector<float>> end;
      std::size_t records_in_block = 0;
      std::string current_chromosome;
      while (it != end && start_pos >= 0)
      {
        if (records_in_block > 0 && it->chromosome() != current_chromosome)
        {
          // TODO: Possibly make this an error case.
          s1r::tree_base::entry e(min, max, (static_cast<std::uint64_t>(start_pos) << 16) | std::uint16_t(records_in_block - 1));
          index_data[current_chromosome].emplace_back(std::move(e));
          min = std::numeric_limits<std::uint32_t>::max();
          max = 0;
        }

        ++records_in_block;
        current_chromosome = it->chromosome();
        min = std::min(min, std::uint32_t(it->locus()));
        max = std::max(max, std::uint32_t(it->locus() + std::max(it->ref().size(), it->alt().size()) - 1));

        std::int64_t end_pos = r.tellg();
        if (start_pos != end_pos) // zstd frame frame boundary
        {
          if (records_in_block > 0x10000) // Max records per block: 64*1024
          {
            assert(!"Too many records in zstd frame to be indexed!");
            return false;
          }

          if (start_pos > 0x0000FFFFFFFFFFFF) // Max file size: 256 TiB
          {
            assert(!"File size to large to be indexed!");
            return false;
          }

          s1r::tree_base::entry e(min, max, (static_cast<std::uint64_t>(start_pos) << 16) | std::uint16_t(records_in_block - 1));
          index_data[it->chromosome()].emplace_back(std::move(e));

          records_in_block = 0;
          min = std::numeric_limits<std::uint32_t>::max();
          max = 0;
        }
        start_pos = end_pos;
        ++it;
      }

      if (start_pos < 0)
      {
        // TODO: handle error.
      }
      else
      {
        ret = s1r::create_file(output_file_path, index_data, 32 - 1);
      }

      return ret;
    }
    //================================================================//
  }
}

