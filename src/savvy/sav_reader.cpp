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
    reader_base::reader_base(const std::string& file_path) :
      subset_size_(0),
      input_stream_(savvy::detail::make_unique<shrinkwrap::zstd::istream>(file_path)),
      file_path_(file_path),
      file_data_format_(fmt::allele)
    {
      parse_header();
      requested_data_format_ = file_data_format_;
    }

    reader_base::reader_base(const std::string& file_path, savvy::fmt data_format) :
      subset_size_(0),
      input_stream_(savvy::detail::make_unique<shrinkwrap::zstd::istream>(file_path)),
      file_path_(file_path),
      file_data_format_(fmt::allele),
      requested_data_format_(data_format)
    {
      parse_header();
    }

    reader_base::reader_base(reader_base&& source) :
      sample_ids_(std::move(source.sample_ids_)),
      subset_map_(std::move(source.subset_map_)),
      subset_size_(source.subset_size_),
      //sbuf_(std::move(source.sbuf_)),
      //input_stream_(&sbuf_),
      input_stream_(std::move(source.input_stream_)),
      file_path_(std::move(source.file_path_)),
      metadata_fields_(std::move(source.metadata_fields_)),
      file_data_format_(source.file_data_format_),
      requested_data_format_(source.requested_data_format_)
    {
    }

    reader_base& reader_base::operator=(reader_base&& source)
    {
      if (&source != this)
      {
        sample_ids_ = std::move(source.sample_ids_);
        subset_map_ = std::move(source.subset_map_);
        subset_size_ = source.subset_size_;
        //sbuf_ = std::move(source.sbuf_);
        //input_stream_->rdbuf(&sbuf_);
        input_stream_ = std::move(source.input_stream_);
        file_path_ = std::move(source.file_path_);
        metadata_fields_ = std::move(source.metadata_fields_);
        file_data_format_ = source.file_data_format_;
        requested_data_format_ = source.requested_data_format_;
      }
      return *this;
    }

    void reader_base::parse_header()
    {
      std::string version_string(7, '\0');
      input_stream_->read(&version_string[0], version_string.size());

      std::string uuid(16, '\0');
      input_stream_->read(&uuid[0], uuid.size());


      std::istreambuf_iterator<char> in_it(*input_stream_);
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
              input_stream_->read(&key[0], key_size);

              std::uint64_t val_size;
              if (varint_decode(in_it, end, val_size) != end)
              {
                ++in_it;
                if (key_size)
                {
                  std::string val;
                  val.resize(val_size);
                  input_stream_->read(&val[0], val_size);

                  if (key == "INFO")
                  {
                    std::string info_field = parse_header_id(val);
                    metadata_fields_.push_back(std::move(info_field));
                  }
                  else if (key == "FORMAT")
                  {
                    std::string format_field = parse_header_id(val);
                    if (format_field == "GT")
                      file_data_format_ = fmt::allele;
//                    else if (format_field == "GP")
//                      file_data_format_ = fmt::genotype_probability;
                    else if (format_field == "HDS")
                      file_data_format_ = fmt::haplotype_dosage;
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
                input_stream_->read(&sample_ids_.back()[0], id_sz);
              }
              --sample_size;
            }

            if (!sample_size)
              return;
          }
        }
      }

      input_stream_->peek();
    }

    std::vector<std::string> reader_base::subset_samples(const std::set<std::string>& subset)
    {
      std::vector<std::string> ret;
      ret.reserve(std::min(subset.size(), sample_ids_.size()));

      subset_map_.clear();
      subset_map_.resize(sample_ids_.size(), std::numeric_limits<std::uint64_t>::max());
      std::uint64_t subset_index = 0;
      for (auto it = sample_ids_.begin(); it != sample_ids_.end(); ++it)
      {
        if (subset.find(*it) != subset.end())
        {
          subset_map_[std::distance(sample_ids_.begin(), it)] = subset_index;
          ret.push_back(*it);
          ++subset_index;
        }
      }

      subset_size_ = subset_index;

      return ret;
    }
    //================================================================//

    //================================================================//
    const std::array<std::string, 0> writer::empty_string_array = {};
    const std::array<std::pair<std::string, std::string>, 0> writer::empty_string_pair_array = {};

    bool writer::create_index(const std::string& input_file_path, std::string output_file_path)
    {
      bool ret = false;

      if (output_file_path.empty())
        output_file_path = input_file_path + ".s1r";

      reader r(input_file_path, fmt::allele); // TODO: make zero if possible.
      std::int64_t start_pos = r.tellg();

      std::uint32_t min = std::numeric_limits<std::uint32_t>::max();
      std::uint32_t max = 0;
      std::map<std::string, std::vector<s1r::tree_base::entry>> index_data;

      site_info variant;
      std::vector<float> variant_data;

      std::size_t records_in_block = 0;
      std::string current_chromosome;
      while (r.read(variant, variant_data) && start_pos >= 0)
      {
        if (records_in_block > 0 && variant.chromosome() != current_chromosome)
        {
          // TODO: Possibly make this an error case.
          s1r::tree_base::entry e(min, max, (static_cast<std::uint64_t>(start_pos) << 16) | std::uint16_t(records_in_block - 1));
          index_data[current_chromosome].emplace_back(std::move(e));
          records_in_block = 0;
          min = std::numeric_limits<std::uint32_t>::max();
          max = 0;
          start_pos = r.tellg();
        }

        ++records_in_block;
        current_chromosome = variant.chromosome();
        min = std::min(min, std::uint32_t(variant.position()));
        max = std::max(max, std::uint32_t(variant.position() + std::max(variant.ref().size(), variant.alt().size()) - 1));

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
          index_data[variant.chromosome()].emplace_back(std::move(e));

          records_in_block = 0;
          min = std::numeric_limits<std::uint32_t>::max();
          max = 0;
          start_pos = end_pos;
        }
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

