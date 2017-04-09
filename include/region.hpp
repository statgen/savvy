
#ifndef LIBVC_REGION_HPP
#define LIBVC_REGION_HPP

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
}

#endif //LIBVC_REGION_HPP
