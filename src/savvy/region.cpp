/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/region.hpp"

namespace savvy
{
  bool region_compare(bounding_point bounding_type, const site_info& var, const region& reg)
  {
    switch (bounding_type)
    {
      case bounding_point::any:   return detail::any_coordinate_within_region::compare(var, reg);
      case bounding_point::all:   return detail::all_coordinates_within_region::compare(var, reg);
      case bounding_point::beg:  return detail::leftmost_coordinate_within_region::compare(var, reg);
      case bounding_point::end: return detail::rightmost_coordinate_within_region::compare(var, reg);
    }
  }
}