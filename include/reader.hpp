#ifndef VC_READER_HPP
#define VC_READER_HPP

#include "haplotype_vector.hpp"
#include "cmf_reader.hpp"
#include "vcf_reader.hpp"
#include "vc.hpp"

#include <string>
#include <memory>
#include <stdexcept>

namespace vc
{

  class reader
  {
  public:
    reader(const std::string& file_path)
    {
      if (vc::detail::has_extension(file_path, ".cmf"))
        cmf_reader_ = std::make_unique<cmf::reader>(file_path);
      else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
        vcf_reader_ = std::make_unique<vcf::reader>(file_path);
    }

    template <typename T>
    reader& operator>>(haplotype_vector<T>& destination)
    {
      if (cmf_reader_)
        (*cmf_reader_) >> destination;
      else if (vcf_reader_)
        (*vcf_reader_) >> destination;
      return *this;
    }

    explicit operator bool() const { return this->good(); }

    bool good() const
    {
      if (cmf_reader_)
        return cmf_reader_->good();
      else if (vcf_reader_)
        return vcf_reader_->good();
      return false;
    }

    bool fail() const
    {
      if (cmf_reader_)
        return cmf_reader_->fail();
      else if (vcf_reader_)
        return vcf_reader_->fail();
      return true;
    }
    bool bad() const
    {
      if (cmf_reader_)
        return cmf_reader_->bad();
      else if (vcf_reader_)
        return vcf_reader_->bad();
      return true;
    }
  private:
    std::unique_ptr<cmf::reader> cmf_reader_;
    std::unique_ptr<vcf::reader> vcf_reader_;
  };



}
#endif //VC_READER_HPP
