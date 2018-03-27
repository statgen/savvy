/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/vcf_reader.hpp"


#include <htslib/vcf.h>

#include <sstream>
#include <assert.h>
#include <limits>

namespace savvy
{
  namespace vcf
  {
    // https://github.com/samtools/htslib/blob/master/tabix.c#L229-L262
    std::vector<std::string> query_chromosomes(const std::string& file_path)
    {
      std::vector<std::string> ret;

      if (::savvy::detail::has_extension(file_path, ".vcf") || ::savvy::detail::has_extension(file_path, ".vcf.gz"))
      {
        tbx_t *tbx = tbx_index_load(file_path.c_str());
        if (tbx)
        {
          const char **seq;
          int nseq{};
          seq = tbx_seqnames(tbx, &nseq);
          ret.resize(nseq);
          for (int i = 0; i < nseq; ++i)
            ret[i] = seq[i];
          tbx_destroy(tbx);
          free(seq);
        }
      }
      else if (::savvy::detail::has_extension(file_path, ".bcf"))
      {
        htsFile *fp = hts_open(file_path.c_str(),"r");
        if (fp)
        {
          bcf_hdr_t *hdr = bcf_hdr_read(fp);
          if (hdr)
          {
            hts_idx_t *idx = bcf_index_load(file_path.c_str());
            if (idx)
            {
              const char **seq;
              int nseq{};
              seq = bcf_index_seqnames(idx, hdr, &nseq);
              ret.resize(nseq);
              for (int i = 0; i < nseq; ++i)
                ret[i] = seq[i];
              hts_idx_destroy(idx);
              free(seq);
            }
            bcf_hdr_destroy(hdr);
          }
          hts_close(fp);
        }
      }

      return ret;
    }
  }
}

//savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a)
//{
//  savvy::vcf::marker::const_iterator ret(a);
//  return (ret += n);
//}
