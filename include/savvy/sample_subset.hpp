/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#if 0
#ifndef LIBSAVVY_SAMPLE_SUBSET_HPP
#define LIBSAVVY_SAMPLE_SUBSET_HPP

#include <vector>
#include <string>

namespace savvy
{
  class sample_subset
  {
  private:
    std::vector<std::string> ids_;
    std::vector<std::size_t> mask_;
  public:
    sample_subset(std::vector<std::string>&& ids, std::vector<std::size_t>&& mask) :
      ids_(std::move(ids)),
      mask_(std::move(mask))
    {
    }

    const std::vector<std::string>& ids() const { return ids_; }
    const std::vector<std::size_t>& mask() const { return mask_; }
  };
}

#endif // LIBSAVVY_SAMPLE_SUBSET_HPP
#endif