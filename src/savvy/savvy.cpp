/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/savvy.hpp"
namespace savvy
{
  namespace detail
  {
    bool has_extension(const std::string& fullString, const std::string& ext)
    {
      if (fullString.length() >= ext.length())
        return (0 == fullString.compare (fullString.length() - ext.length(), ext.length(), ext));
      else
        return false;
    }
  }
}