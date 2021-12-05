/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_LOGGING_HPP
#define LIBSAVVY_LOGGING_HPP

#include <unordered_set>
#include <cstdio>

template <typename T = void>
class logging
{
private:
 static std::unordered_set<std::string> distinct_messages_;
public:
  template<typename... A>
  static void cerr_once(const std::string& s, A ...args)
  {
    std::string buf(512, '\0');
    int sz = std::sprintf(&buf[0], s.c_str(), args...);
    if (sz < 0)
      return std::cerr << "Warning: Log failed\n", void();
    else if (sz == (int)buf.size())
      std::cerr << "Warning: log message too long\n";

    buf.resize(sz);
    if (distinct_messages_.insert(buf).second)
      std::cerr.write(buf.data(), buf.size());
  }
};

template <typename T>
std::unordered_set<std::string> logging<T>::distinct_messages_;

#endif // LIBSAVVY_LOGGING_HPP