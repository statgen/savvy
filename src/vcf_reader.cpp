
#include "vcf_reader.hpp"

#include <assert.h>

namespace vc
{
  namespace vcf
  {
    marker::marker() :
      hts_rec_(bcf_init1())
    {

    }

    marker::~marker()
    {
      try
      {
        if (hts_rec_)
          bcf_destroy1(hts_rec_);
      }
      catch (...)
      {
        assert(!"bcf_destroy1 threw exception!");
      }
    }

    bool marker::read_marker(marker& destination, htsFile* hts_file, bcf_hdr_t* hts_hdr)
    {
      return (bcf_read(hts_file, hts_hdr, destination.hts_rec_) >= 0);
    }

    reader::reader(const std::string& file_path) :
      hts_file_(bcf_open(file_path.c_str(), "r")), hts_hdr_(nullptr)
    {
      if (hts_file_)
      {
        hts_hdr_ = vcf_hdr_read(hts_file_);
      }
    }

    reader::~reader()
    {
      try
      {
        if (hts_hdr_)
          bcf_hdr_destroy(hts_hdr_);
        if (hts_file_)
          bcf_close(hts_file_);
      }
      catch (...)
      {
        assert(!"bcf_hdr_destroy or bcf_close threw exception!");
      }
    }

    bool reader::read_next_marker(marker& destination)
    {
      return marker::read_marker(destination, hts_file_, hts_hdr_);
    }
  }
}
