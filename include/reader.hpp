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
  template <typename CMFType, typename VCFType>
  class reader_base
  {
  public:
    operator bool() const
    {
      return this->good();
    }

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
  protected:
    template <typename T>
    void read_vector(haplotype_vector<T>& destination)
    {
      if (cmf_reader_)
        (*cmf_reader_) >> destination;
      else if (vcf_reader_)
        (*vcf_reader_) >> destination;
    }
  protected:
    std::unique_ptr<CMFType> cmf_reader_;
    std::unique_ptr<VCFType> vcf_reader_;
  };

  class reader : public reader_base<cmf::reader, vcf::reader>
  {
  public:
    reader(const std::string& file_path);
    ~reader() {}

    template <typename T>
    reader& operator>>(haplotype_vector<T>& destination);
  };

  class region_reader : public reader_base<cmf::region_reader, vcf::region_reader>
  {
  public:
    region_reader(const std::string& file_path, region reg);

    template <typename T>
    region_reader& operator>>(haplotype_vector<T>& destination);
  private:

  };

  template <typename T>
  reader& reader::operator>>(haplotype_vector<T>& destination)
  {
    read_vector(destination);
    return *this;
  }

  template <typename T>
  region_reader& region_reader::operator>>(haplotype_vector<T>& destination)
  {
    read_vector(destination);
    return *this;
  }
}
#endif //VC_READER_HPP
