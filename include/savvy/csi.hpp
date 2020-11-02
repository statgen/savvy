/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_CSI_HPP
#define LIBSAVVY_CSI_HPP

#include <shrinkwrap/gz.hpp>

#include <array>
#include <unordered_map>
#include <vector>
#include <cstdint>

class csi_index
{
public:
  bool good() const { return fs_.good(); }

  csi_index(const std::string& file_path)
    : fs_(file_path)
  {
    std::array<char, 4> magic;
    fs_.read(magic.data(), magic.size());

    fs_.read((char*)&min_shift_, sizeof(min_shift_));
    fs_.read((char*)&depth_, sizeof(depth_));

    std::uint32_t aux_sz{};
    std::vector<std::uint8_t> aux_data;
    if (fs_.read((char*)&aux_sz, sizeof(aux_sz)))
    {
      aux_data.resize(aux_sz);
      fs_.read((char*)aux_data.data(), aux_sz);

      std::size_t off = 28;
      while (off < aux_sz)
      {
        const char* p = (char*)aux_data.data() + off;
        aux_contigs_.emplace_back(std::string(p));
        off += aux_contigs_.back().size() + 1;
      }

      std::uint32_t n_indices{};
      if (fs_.read((char*)&n_indices, sizeof(n_indices)))
      {
        indices_.resize(n_indices);
        for (std::size_t i = 0; i < n_indices; ++i)
        {
          std::uint32_t n_bins{};
          if (fs_.read((char*)&n_bins, sizeof(n_bins)))
          {
            indices_[i].reserve(n_bins);
            for (std::size_t b = 0; b < n_bins; ++b)
            {
              std::uint32_t bin_id{}, n_chunks{};
              std::uint64_t voff{};

              fs_.read((char*)&bin_id, sizeof(bin_id));
              fs_.read((char*)&voff, sizeof(voff));

              if (fs_.read((char*)&n_chunks, sizeof(n_chunks)))
              {
                auto& bin = indices_[i][bin_id];
                bin.loff = voff;
                bin.chunks.resize(n_chunks);
                for (std::size_t c = 0; c < n_chunks; ++c)
                {
                  fs_.read((char*)&bin.chunks[c].first, 8);
                  fs_.read((char*)&bin.chunks[c].second, 8);
                }
              }
            }
          }
        }
      }
    }

    if (!fs_)
    {
      std::cerr << "Error: malformed csi index" << std::endl;
    }

  }

  /* calculate maximum bin number -- valid bin numbers range within [0,bin_limit) */
  int bin_limit()
  {
    return ((1 << (depth_+1)*3) - 1) / 7;
  }

  /* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
  std::vector<int> reg2bins(std::int64_t beg, std::int64_t end)
  {
    std::vector<int> bins;
    int l, t, s = min_shift_ + depth_*3;
    for (--end, l = t = 0; l <= depth_; s -= 3, t += 1<<l*3, ++l) {
      int b = t + (beg>>s), e = t + (end>>s), i;
      for (i = b; i <= e; ++i) bins.emplace_back(i);
    }
    return bins;
  }

  class bin_t
  {
  public:
    std::uint64_t loff = 0;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> chunks;
  };

  std::vector<std::pair<std::uint64_t, std::uint64_t>> query_intervals(const std::string& contig, const std::unordered_map<std::string, std::uint32_t>& contig_to_id, std::int64_t beg, std::int64_t end)
  {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> ret;
    auto bin_ids = reg2bins(beg, end);
    ret.reserve(bin_ids.size());

    auto contig_id = std::uint32_t(-1);
    if (aux_contigs_.size())
    {
      for (std::size_t i = 0; i < aux_contigs_.size(); ++i)
      {
        if (aux_contigs_[i] == contig)
        {
          contig_id = i;
          break;
        }
      }
    }
    else
    {
      auto it = contig_to_id.find(contig);
      if (it != contig_to_id.end())
        contig_id = it->second;
    }

    if (contig_id < indices_.size())
    {
      ret.reserve(bin_ids.size());
      for (auto bin_id_it = bin_ids.begin(); bin_id_it != bin_ids.end(); ++bin_id_it)
      {
        auto bin_it = indices_[contig_id].find(*bin_id_it);
        if (bin_it != indices_[contig_id].end())
        {
          assert(bin_it->first < bin_limit());
          for (auto chunk_it = bin_it->second.chunks.begin(); chunk_it != bin_it->second.chunks.end(); ++chunk_it)
          {
            if (ret.empty())
            {
              ret.emplace_back(chunk_it->first, chunk_it->second);
            }
            else
            {
              // reversed insertion sort
              auto it = ret.end();

              do
              {
                --it;
                if (it->first <= chunk_it->first)
                {
                  ret.emplace(it + 1, chunk_it->first, chunk_it->second);
                  break;
                }

                if (it == ret.begin())
                  it = ret.emplace(ret.begin(), chunk_it->first, chunk_it->second);

              } while (it != ret.begin());
            }
          }
        }
      }

      // merge intervals
      if (ret.size() > 0)
      {
        auto cur = ret.begin();
        auto next = cur + 1;

        while(next != ret.end())
        {
          if(next->second < cur->second) // next fits inside cur
          {
            next = ret.erase(next);
            cur = next - 1;
          }
          else if((next->first >> 16) <= (cur->second >> 16)) // next overlaps cur
          {
            cur->second = next->second;
            next = ret.erase(next);
            cur = next - 1;
          }
          else
          {
            cur = next;
            ++next;
          }
        }
      }
    }

    return ret;
  }


private:
  shrinkwrap::bgzf::istream fs_;
  std::int32_t min_shift_ = 0;
  std::int32_t depth_ = 0;
  std::vector<std::string> aux_contigs_;
  std::vector<std::unordered_map<std::uint32_t, bin_t>> indices_;
};

#endif // LIBSAVVY_CSI_HPP