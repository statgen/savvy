#ifndef LIBSAVVY_UTILITY_HPP
#define LIBSAVVY_UTILITY_HPP

#include <string>
#include <memory>

namespace savvy
{
//  template<typename T, typename... Args>
//  std::unique_ptr<T> make_unique(Args&&... args)
//  {
//    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
//  }

  std::string parse_header_id(std::string header_value);
}

#endif // LIBSAVVY_UTILITY_HPP