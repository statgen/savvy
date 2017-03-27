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
  template <typename VecType>
  class reader_impl_base
  {
  public:
    virtual operator bool() const = 0;
    virtual bool good() const = 0;
    virtual bool fail() const = 0;
    virtual bool bad() const = 0;
    virtual void read(haplotype_vector<VecType>& destination) = 0;
  };

  template <typename VecType>
  class cmf_reader_impl : public reader_impl_base<VecType>
  {
  public:
    cmf_reader_impl(const std::string& file_path) : r_(file_path) {}
    operator bool() const { return this->good(); }
    bool good() const { return this->r_.good(); }
    bool fail() const { return this->r_.fail(); }
    bool bad() const { return this->r_.bad(); }
    void read(haplotype_vector<VecType>& destination)
    {
      r_.read(destination);
    }
  private:
    cmf::reader r_;
  };

  template <typename VecType>
  class cmf_region_reader_impl : public reader_impl_base<VecType>
  {
  public:
    cmf_region_reader_impl(const std::string& file_path, region reg) : r_(file_path, reg) {}
    operator bool() const { return this->good(); }
    bool good() const { return this->r_.good(); }
    bool fail() const { return this->r_.fail(); }
    bool bad() const { return this->r_.bad(); }
    void read(haplotype_vector<VecType>& destination)
    {
      r_.read(destination);
    }
  private:
    cmf::region_reader r_;
  };

  template <typename VecType>
  class vcf_reader_impl : public reader_impl_base<VecType>
  {
  public:
    vcf_reader_impl(const std::string& file_path) : r_(file_path) {}
    operator bool() const { return this->good(); }
    bool good() const { return this->r_.good(); }
    bool fail() const { return this->r_.fail(); }
    bool bad() const { return this->r_.bad(); }
    void read(haplotype_vector<VecType>& destination)
    {
      r_.read(destination);
    }
  private:
    vcf::reader r_;
  };

  template <typename VecType>
  class vcf_region_reader_impl : public reader_impl_base<VecType>
  {
  public:
    vcf_region_reader_impl(const std::string& file_path, region reg) : r_(file_path, reg) {}
    operator bool() const { return this->good(); }
    bool good() const { return this->r_.good(); }
    bool fail() const { return this->r_.fail(); }
    bool bad() const { return this->r_.bad(); }
    void read(haplotype_vector<VecType>& destination)
    {
      r_.read(destination);
    }
  private:
    vcf::region_reader r_;
  };


  template <typename VecType>
  class reader_base
  {
  public:
    operator bool() const
    {
      return this->good();
    }

    bool good() const
    {
      if (impl_)
        return impl_->good();
      return false;
    }

    bool fail() const
    {
      if (impl_)
        return impl_->fail();
      return true;
    }

    bool bad() const
    {
      if (impl_)
        return impl_->bad();
      return true;
    }

    void read(haplotype_vector<VecType>& destination)
    {
      if (impl_)
        impl_->read(destination);
    }
  protected:
    std::unique_ptr<reader_impl_base<VecType>> impl_;
  };

  template <typename VecType>
  class basic_reader : public reader_base<VecType>
  {
  public:
    basic_reader(const std::string& file_path);
    ~basic_reader() {}

    basic_reader& operator>>(haplotype_vector<VecType>& destination);
  };

  template <typename VecType>
  class basic_region_reader : public reader_base<VecType>
  {
  public:
    basic_region_reader(const std::string& file_path, region reg);

    basic_region_reader& operator>>(haplotype_vector<VecType>& destination);
  };

  template <typename T>
  basic_reader<T>::basic_reader(const std::string& file_path)
  {
    if (vc::detail::has_extension(file_path, ".cmf"))
      reader_base<T>::impl_ = std::make_unique<cmf_reader_impl<T>>(file_path);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      reader_base<T>::impl_ = std::make_unique<vcf_reader_impl<T>>(file_path);
  }

  template <typename T>
  basic_reader<T>& basic_reader<T>::operator>>(haplotype_vector<T>& destination)
  {
    read(destination);
    return *this;
  }

  template <typename T>
  basic_region_reader<T>::basic_region_reader(const std::string& file_path, region reg)
  {
    if (vc::detail::has_extension(file_path, ".cmf"))
      reader_base<T>::impl_ = std::make_unique<cmf_region_reader_impl<T>>(file_path, reg);
    else if (detail::has_extension(file_path, ".vcf") || detail::has_extension(file_path, ".vcf.gz") || detail::has_extension(file_path, ".bcf"))
      reader_base<T>::impl_ = std::make_unique<vcf_region_reader_impl<T>>(file_path, reg);
  }

  template <typename T>
  basic_region_reader<T>& basic_region_reader<T>::operator>>(haplotype_vector<T>& destination)
  {
    read(destination);
    return *this;
  }

  typedef basic_reader<std::vector<float>> reader;
  typedef basic_region_reader<std::vector<float>> region_reader;

  using variant_iterator =  basic_variant_iterator<reader_base<std::vector<float>>, std::vector<float>>;
}
#endif //VC_READER_HPP
