/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/region.hpp"

namespace savvy
{
  bool region_compare(coord_bound bounding_type, const site_info& var, const region& reg)
  {
    switch (bounding_type)
    {
      case coord_bound::any:   return detail::any_coordinate_within_region::compare(var, reg);
      case coord_bound::all:   return detail::all_coordinates_within_region::compare(var, reg);
      case coord_bound::left:  return detail::leftmost_coordinate_within_region::compare(var, reg);
      case coord_bound::right: return detail::rightmost_coordinate_within_region::compare(var, reg);
    }
  }
}