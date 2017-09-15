
#ifndef LIBSAVVY_REGION_HPP
#define LIBSAVVY_REGION_HPP

#include "site_info.hpp"

#include <cstdint>
#include <string>
#include <limits>

namespace savvy
{
  class region
  {
  public:
    region(const std::string& chromosome, std::uint64_t from = 1, std::uint64_t to = std::numeric_limits<std::uint64_t>::max()) :
      chromosome_(chromosome),
      from_(from),
      to_(to)
    {
    }
    const std::string& chromosome() const { return chromosome_; }
    std::uint64_t from() const { return from_; }
    std::uint64_t to() const { return to_; }
  private:
    std::string chromosome_;
    std::uint64_t from_;
    std::uint64_t to_;
  };

  struct any_coordinate_within_region
  {
    bool operator()(const site_info& var, const region& reg)
    {
      return (var.locus() <= reg.to() && (var.locus() + std::max(var.ref().size(), var.alt().size()) - 1) >= reg.from() && var.chromosome() == reg.chromosome());
    }
  };

  struct all_coordinates_within_region
  {
    bool operator()(const site_info& var, const region& reg)
    {
      return (var.locus() >= reg.from() && (var.locus() + std::max(var.ref().size(), var.alt().size()) - 1) <= reg.to() && var.chromosome() == reg.chromosome());
    }
  };

  struct leftmost_coordinate_within_region
  {
    bool operator()(const site_info& var, const region& reg)
    {
      return (var.locus() >= reg.from() && var.locus() <= reg.to() && var.chromosome() == reg.chromosome());
    }
  };

  struct rightmost_coordinate_within_region
  {
    bool operator()(const site_info& var, const region& reg)
    {
      std::uint64_t right = (var.locus() + std::max(var.ref().size(), var.alt().size()) - 1);
      return (right >= reg.from() && right <= reg.to() && var.chromosome() == reg.chromosome());
    }
  };
}

#endif //LIBSAVVY_REGION_HPP
