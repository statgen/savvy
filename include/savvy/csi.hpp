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

namespace  savvy
{
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
    std::vector<int> reg2bins_old(std::int64_t beg, std::int64_t end)
    {
      std::vector<int> bins;
      int l, t, s = min_shift_ + depth_*3;
      for (--end, l = t = 0; l <= depth_; s -= 3, t += 1<<l*3, ++l) {
        int b = t + (beg>>s), e = t + (end>>s), i;
        for (i = b; i <= e; ++i) bins.emplace_back(i);
      }
      return bins;
    }

    std::vector<int> reg2bins(std::int64_t beg, std::int64_t end)
    {
      std::vector<int> bins;
      int l, t, s = min_shift_ + depth_ * 3;
      if (beg >= end) return bins;
      if (end >= 1LL<<s) end = 1LL<<s;
      for (--end, l = 0, t = 0; l <= depth_; s -= 3, t += 1<<(l*3), ++l) {
        std::int64_t b, e;
        //int n, i;
        b = t + (beg>>s); e = t + (end>>s); //n = e - b + 1;
        for (int i = b; i <= e; ++i) bins.emplace_back(i);
      }
      return bins;
    }

    class bin_t
    {
    public:
      std::uint64_t loff = 0;
      std::vector<std::pair<std::uint64_t, std::uint64_t>> chunks;
    };

    static int bin_first(int l) { return (((1<<(((l)<<1) + (l))) - 1) / 7); }
    static int bin_parent(int l) { return (((l) - 1) >> 3); }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> query_intervals(const std::string& contig, const std::unordered_map<std::string, std::uint32_t>& contig_to_id, std::int64_t beg, std::int64_t end)
    {
      std::vector<std::pair<std::uint64_t, std::uint64_t>> ret;
      auto bin_ids = reg2bins(beg, end);

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
#if 1
        //===============================================//
        // Experimental
        // Adapted from https://github.com/samtools/htslib/blob/a7a90fe913f8a466f32f6e284cf46653944acd6f/hts.c#L2602-L2684

        // compute min_off
        int bin = bin_first(depth_) + (beg>>min_shift_);
        typename decltype(indices_)::value_type::iterator bin_it;
        do {
          int first;
          bin_it = indices_[contig_id].find(bin); //k = kh_get(bin, bidx, bin);
          if (bin_it != indices_[contig_id].end())
            break;
          first = (bin_parent(bin)<<3) + 1;
          if (bin > first) --bin;
          else bin = bin_parent(bin);
        } while (bin);
        if (bin == 0)
          bin_it = indices_[contig_id].find(bin);
        std::uint64_t min_off = bin_it != indices_[contig_id].end() ? bin_it->second.loff : 0;
//        if (idx->lidx[tid].offset
//          && beg>>idx->min_shift < idx->lidx[tid].n
//          && min_off < idx->lidx[tid].offset[beg>>idx->min_shift])
//          min_off = idx->lidx[tid].offset[beg>>idx->min_shift];

        // compute max_off: a virtual offset from a bin to the right of end
        std::uint64_t max_off = 0;
        bin = bin_first(depth_) + ((end-1) >> min_shift_) + 1;
        if (bin >= bin_limit()) bin = 0;
        while (1)
        {
          // search for an extant bin by moving right, but moving up to the
          // parent whenever we get to a first child (which also covers falling
          // off the RHS, which wraps around and immediately goes up to bin 0)
          while (bin % 8 == 1)
            bin = bin_parent(bin);

          if (bin == 0)
          {
            max_off = (uint64_t)-1;
            break;
          }
          bin_it = indices_[contig_id].find(bin);
          if (bin_it != indices_[contig_id].end() && bin_it->second.chunks.size() > 0)
          {
            max_off = bin_it->second.chunks[0].first;
            break;
          }
          bin++;
        }


        std::size_t n_off = 0;
        for (std::size_t i = 0; i < bin_ids.size(); ++i)
        {
          auto bin_it = indices_[contig_id].find(bin_ids[i]);
          if (bin_it != indices_[contig_id].end()) //if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx))
            n_off += bin_it->second.chunks.size(); //kh_value(bidx, k).n;
        }
        if (n_off == 0)
        {
          // No overlapping bins means the iterator has already finished.
          //iter->finished = 1;
          return ret; //iter;
        }
        //off = calloc(n_off, sizeof(*off));
        for (std::size_t i = n_off = 0; i < bin_ids.size(); ++i)
        {
          auto bin_it = indices_[contig_id].find(bin_ids[i]);
          if (bin_it != indices_[contig_id].end()) //if ((k = kh_get(bin, bidx, iter->bins.a[i])) != kh_end(bidx)) {
          {
            //bins_t *p = &kh_value(bidx, k);
            auto& list = bin_it->second.chunks;
            for (std::size_t j = 0; j < list.size(); ++j)
            {
              if (list[j].second > min_off && list[j].first < max_off)
              {
//                off[n_off].u = min_off > list[j].u
//                               ? min_off : list[j].u;
//                off[n_off].v = max_off < list[j].v
//                               ? max_off : p->list[j].v;
                ret.emplace_back(
                  min_off > list[j].first ? min_off : list[j].first,
                  max_off < list[j].second ? max_off : list[j].second);
                // hts_pair64_max_t::max is now used to link
                // file offsets to region list entries.
                // The iterator can use this to decide if it
                // can skip some file regions.
                //off[n_off].max = ((uint64_t) tid << 32) | j;
                n_off++;
              }
            }
          }
        }

        if (n_off == 0)
        {
          assert(ret.empty());
          return ret;
        }
        //ks_introsort(_off_max, n_off, off);
        std::sort(ret.begin(), ret.end(), [](const std::pair<std::uint64_t, std::uint64_t>& l, const std::pair<std::uint64_t, std::uint64_t>& r)
        {
          return l.first < r.first;
        });
        // resolve completely contained adjacent blocks
        assert(n_off <= ret.size());
        std::size_t l = 0;
        for (std::size_t i = 1; i < n_off; ++i)
        {
          if (ret[l].second < ret[i].second)
            ret[++l] = ret[i];
        }
        n_off = l + 1;
        // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
        for (std::size_t i = 1; i < n_off; ++i)
          if (ret[i-1].second >= ret[i].first) ret[i-1].second = ret[i].first;
        // merge adjacent blocks
        l = 0;
        for (std::size_t i = 1; i < n_off; ++i)
        {
          if ((ret[l].second >> 16) == (ret[i].first >> 16))
            ret[l].second = ret[i].second;
          else
            ret[++l] = ret[i];
        }
        n_off = l + 1;
        assert(n_off <= ret.size());
        ret.resize(n_off);
        //===============================================//
#else
        for (auto bin_id_it = bin_ids.begin(); bin_id_it != bin_ids.end(); ++bin_id_it)
        {
          auto bin_it = indices_[contig_id].find(*bin_id_it);
          if (bin_it != indices_[contig_id].end())
          {
            assert(bin_it->first < std::size_t(bin_limit()));
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
#endif
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
}

#endif // LIBSAVVY_CSI_HPP