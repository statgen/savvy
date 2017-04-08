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
    bool read(haplotype_vector<T>& destination)
    {
      if (cmf_reader_)
        cmf_reader_->read(destination);
      else if (vcf_reader_)
        vcf_reader_->read(destination);

      return good();
    }
  protected:
    std::unique_ptr<sav::reader_base> cmf_reader_;
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

  class indexed_reader : public reader_base
  {
  public:
    indexed_reader(const std::string& file_path, const region& reg);
    void reset_region(const region& reg);

    template <typename T>
    indexed_reader& operator>>(haplotype_vector<T>& destination);
    template <typename T, typename Pred>
    indexed_reader& read_if(haplotype_vector<T>& destination, Pred fn, const typename T::value_type missing_value = std::numeric_limits<typename T::value_type>::quiet_NaN(), const typename T::value_type alt_value = 1, const typename T::value_type ref_value = 0);
  private:

  };

  template <typename T>
  reader& reader::operator>>(haplotype_vector<T>& destination)
  {
    read(destination);
    return *this;
  }

  template <typename T>
  indexed_reader& indexed_reader::operator>>(haplotype_vector<T>& destination)
  {
    read(destination);
    return *this;
  }

  template <typename T, typename Pred>
  indexed_reader& indexed_reader::read_if(haplotype_vector<T>& destination, Pred fn, const typename T::value_type missing_value, const typename T::value_type alt_value, const typename T::value_type ref_value)
  {
    if (cmf_reader_)
      dynamic_cast<sav::indexed_reader*>(cmf_reader_.get())->read_if(destination, fn, missing_value, alt_value, ref_value);
    else if (vcf_reader_)
      dynamic_cast<vcf::indexed_reader*>(vcf_reader_.get())->read_if(destination, fn, missing_value, alt_value, ref_value);

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