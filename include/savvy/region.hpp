/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_REGION_HPP
#define LIBSAVVY_REGION_HPP

#include "site_info.hpp"

#include <cstdint>
#include <string>
#include <limits>
#include <vector>
#include <unordered_map>

namespace savvy
{
  /// Bounding point enum.
  enum class bounding_point : std::uint8_t
  {
    any = 0, ///< Variants for which any nucleotide in reference or alternate alleles are within query region
    all, ///< Variants for which all nucleotides in reference or alternate alleles are within query region
    beg, ///< Variants for which begin position are within query region
    end ///< Variants for which end position (i.e., position + max[size_of_ref, size_of_alts]) are within query region
  };

  class query_bounds
  {
  public:
    /**
     * Merges container of overlapping regions.
     * @tparam Iter Iterator type
     * @param beg Begin iterator of container
     * @param end End iterator of container
     * @return Vector of merged regions
     */
    template <typename Iter>
    static std::vector<query_bounds> merge(Iter beg, Iter end);

    query_bounds(const std::string& chromosome, std::uint64_t from = 1, std::uint64_t to = std::numeric_limits<std::uint64_t>::max()) :
      chromosome_(chromosome),
      from_(from),
      to_(to)
    {
    }

    /**
     * Gets chromosome for query.
     * @return Chromosome string
     */
    const std::string& chromosome() const { return chromosome_; }

    /**
     * Gets start position of query.
     * @return Start position
     */
    std::uint64_t from() const { return from_; }

    /**
     * Gets end position of query. Genomic regions are 1-based and inclusive end. Slice bounds are 0-based and exclusive end.
     * @return End position
     */
    std::uint64_t to() const { return to_; }
  private:
    std::string chromosome_;
    std::uint64_t from_;
    std::uint64_t to_;
  };

  class genomic_region : public query_bounds
  {
  public:
    genomic_region(query_bounds&& src) : query_bounds(std::move(src)) {}
    genomic_region(const std::string& chromosome, std::uint64_t from = 0, std::uint64_t to = std::numeric_limits<std::uint64_t>::max()) :
      query_bounds(chromosome, from, to)
    {
    }
  };

  typedef genomic_region region;

  class slice_bounds : public query_bounds
  {
  public:
    slice_bounds(query_bounds&& src) : query_bounds(std::move(src)) {}
    slice_bounds(std::uint64_t from, std::uint64_t to = std::numeric_limits<std::uint64_t>::max(), const std::string& chromosome = "") :
      query_bounds(chromosome, from, to)
    {
    }
  };

  template <typename Iter>
  std::vector<query_bounds> query_bounds::merge(Iter beg, Iter end)
  {
    std::unordered_map<std::string, std::size_t> ret_index;
    std::vector<query_bounds> ret;

    for (auto it = beg; it != end; ++it)
    {
      auto insert_res = ret_index.insert(std::make_pair(it->chromosome(), ret.size()));
      if (insert_res.second)
      {
        ret.emplace_back(*it);
      }
      else
      {
        std::uint64_t from = std::min(ret[insert_res.first->second].from(), it->from());
        std::uint64_t to = std::max(ret[insert_res.first->second].to(), it->to());
        ret[insert_res.first->second] = genomic_region(ret[insert_res.first->second].chromosome(), from, to);
      }
    }

    return ret;
  }

  //namespace v2
  //{
    namespace detail
    {
      struct any_coordinate_within_region
      {
        static bool compare(const site_info& var, const genomic_region& reg)
        {
          std::uint32_t end_val;
          if (!var.get_info("END", reinterpret_cast<std::int32_t&>(end_val)))
          {
            std::uint32_t max_allele_size(var.ref().size());
            for (auto it = var.alts().begin(); it != var.alts().end(); ++it)
              max_allele_size = std::max(max_allele_size, std::uint32_t(it->size()));
            end_val = var.pos() + std::max(std::uint32_t(var.ref().size()), max_allele_size) - 1;
          }
          return (var.pos() <= reg.to() && end_val >= reg.from() && (var.chrom() == reg.chromosome() || reg.chromosome().empty()));
        }
      };

      struct all_coordinates_within_region
      {
        static bool compare(const site_info& var, const genomic_region& reg)
        {
          std::uint32_t end_val;
          if (!var.get_info("END", reinterpret_cast<std::int32_t&>(end_val)))
          {
            std::uint32_t max_allele_size(var.ref().size());
            for (auto it = var.alts().begin(); it != var.alts().end(); ++it)
              max_allele_size = std::max(max_allele_size, std::uint32_t(it->size()));
            end_val = var.pos() + std::max(std::uint32_t(var.ref().size()), max_allele_size) - 1;
          }
          return (var.pos() >= reg.from() && end_val <= reg.to() && (var.chrom() == reg.chromosome() || reg.chromosome().empty()));
        }
      };

      struct leftmost_coordinate_within_region
      {
        static bool compare(const site_info& var, const genomic_region& reg)
        {
          return (var.pos() >= reg.from() && var.pos() <= reg.to() && (var.chrom() == reg.chromosome() || reg.chromosome().empty()));
        }
      };

      struct rightmost_coordinate_within_region
      {
        static bool compare(const site_info& var, const genomic_region& reg)
        {
          std::uint32_t end_val;
          if (!var.get_info("END", reinterpret_cast<std::int32_t&>(end_val)))
          {
            std::uint32_t max_allele_size(var.ref().size());
            for (auto it = var.alts().begin(); it != var.alts().end(); ++it)
              max_allele_size = std::max(max_allele_size, std::uint32_t(it->size()));
            end_val = (var.pos() + std::max(std::uint32_t(var.ref().size()), max_allele_size) - 1);
          }
          return (end_val >= reg.from() && end_val <= reg.to() && (var.chrom() == reg.chromosome() || reg.chromosome().empty()));
        }
      };
    }

    inline bool region_compare(bounding_point bounding_type, const site_info& var, const genomic_region& reg)
    {
      switch (bounding_type)
      {
      case bounding_point::any:
        return detail::any_coordinate_within_region::compare(var, reg);
      case bounding_point::all:
        return detail::all_coordinates_within_region::compare(var, reg);
      case bounding_point::beg:
        return detail::leftmost_coordinate_within_region::compare(var, reg);
      case bounding_point::end:
        return detail::rightmost_coordinate_within_region::compare(var, reg);
      default:
        return false;
      }
    }
  //}

#if 0
  namespace detail
  {
    struct any_coordinate_within_region
    {
      static bool compare(const site_info& var, const genomic_region& reg)
      {
        return (var.position() <= reg.to() && (var.position() + std::max(var.ref().size(), var.alt().size()) - 1) >= reg.from() && (var.chromosome() == reg.chromosome() || reg.chromosome().empty()));
      }
    };

    struct all_coordinates_within_region
    {
      static bool compare(const site_info& var, const genomic_region& reg)
      {
        return (var.position() >= reg.from() && (var.position() + std::max(var.ref().size(), var.alt().size()) - 1) <= reg.to() && (var.chromosome() == reg.chromosome() || reg.chromosome().empty()));
      }
    };

    struct leftmost_coordinate_within_region
    {
      static bool compare(const site_info& var, const genomic_region& reg)
      {
        return (var.position() >= reg.from() && var.position() <= reg.to() && (var.chromosome() == reg.chromosome() || reg.chromosome().empty()));
      }
    };

    struct rightmost_coordinate_within_region
    {
      static bool compare(const site_info& var, const genomic_region& reg)
      {
        std::uint64_t right = (var.position() + std::max(var.ref().size(), var.alt().size()) - 1);
        return (right >= reg.from() && right <= reg.to() && (var.chromosome() == reg.chromosome() || reg.chromosome().empty()));
      }
    };
  }

  bool region_compare(bounding_point bounding_type, const site_info& var, const genomic_region& reg);
#endif
}

#endif //LIBSAVVY_REGION_HPP
