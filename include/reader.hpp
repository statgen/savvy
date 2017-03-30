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

    template <typename T>
    void read(haplotype_vector<T>& destination)
    {
      if (cmf_reader_)
        cmf_reader_->read(destination);
      else if (vcf_reader_)
        vcf_reader_->read(destination);
    }
  protected:
    std::unique_ptr<cmf::reader_base> cmf_reader_;
    std::unique_ptr<vcf::reader_base> vcf_reader_;
  };

  class reader : public reader_base
  {
  public:
    reader(const std::string& file_path);
    ~reader() {}

    template <typename T>
    reader& operator>>(haplotype_vector<T>& destination);
  };

  class region_reader : public reader_base
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
    read(destination);
    return *this;
  }

  template <typename T>
  region_reader& region_reader::operator>>(haplotype_vector<T>& destination)
  {
    read(destination);
    return *this;
  }

  template <typename VecType>
  using variant_iterator =  basic_variant_iterator<reader_base, VecType>;

  template <typename ValType>
  using dense_variant_iterator =  basic_variant_iterator<reader_base, std::vector<ValType>>;
  template <typename ValType>
  using sparse_variant_iterator =  basic_variant_iterator<reader_base, compressed_vector<ValType>>;
}

#endif //VC_READER_HPP