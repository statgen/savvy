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

      reader<1> r(input_file_path, fmt::allele); // TODO: make zero if possible.
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
          min = std::numeric_limits<std::uint32_t>::max();
          max = 0;
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
        }
        start_pos = end_pos;
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

