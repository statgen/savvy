/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LIBSAVVY_FILE_HPP
#define LIBSAVVY_FILE_HPP

#include "dictionary.hpp"

namespace savvy
{
  class file
  {
  protected:
    dictionary dict_;
  public:
    enum class format
    {
      sav1 = 1,
      sav2,
      bcf,
      vcf
    };

    const ::savvy::dictionary& dictionary() const { return dict_; }
    virtual ~file() {}
  };
}

#endif // LIBSAVVY_FILE_HPP