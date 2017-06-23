#ifndef VC_READER_HPP
#define VC_READER_HPP

#include "allele_vector.hpp"
#include "sav_reader.hpp"
#include "vcf_reader.hpp"
#include "savvy.hpp"

#include <string>
#include <memory>
#include <stdexcept>

namespace savvy
{
  class reader_base
  {
  public:
    class sample_iterator
    {
    public:
      typedef sample_iterator self_type;
      typedef std::ptrdiff_t difference_type;
      typedef std::string value_type;
      typedef const value_type& reference;
      typedef const value_type* pointer;
      typedef std::bidirectional_iterator_tag iterator_category;
    public:
      sample_iterator() : cstring_itr_(nullptr), stdstring_itr_{} {}
      sample_iterator(const char** cstring_itr) : cstring_itr_(cstring_itr), stdstring_itr_{} {}
      sample_iterator(std::vector<std::string>::const_iterator stdstring_itr) : cstring_itr_(nullptr), stdstring_itr_(stdstring_itr) {}

      self_type& operator+=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ += n;
        else
          stdstring_itr_ += n;
        return *this;
      }
      self_type operator+(difference_type n) const { self_type ret(*this); return (ret += n); }
      self_type& operator-=(difference_type n)
      {
        if (cstring_itr_)
          cstring_itr_ -= n;
        else
          stdstring_itr_ -= n;
        return *this;
      }
      self_type operator-(difference_type n) const { self_type ret(*this); return (ret -= n); }
      difference_type operator-(const self_type& b) const
      {
        if (cstring_itr_)
          return cstring_itr_ - b.cstring_itr_;
        return stdstring_itr_ - b.stdstring_itr_;
      }

      self_type& operator--()
      {
        if (cstring_itr_)
          --cstring_itr_;
        else
          --stdstring_itr_;
        return *this;
      }
      self_type operator--(int) { self_type r = *this; --(*this); return r; }
      self_type& operator++()
      {
        if (cstring_itr_)
          ++cstring_itr_;
        else
          ++stdstring_itr_;
        return *this;
      }
      self_type operator++(int) { self_type r = *this; ++(*this); return r; }
      reference operator*()
      {
        if (cstring_itr_)
        {
          tmp_ = *cstring_itr_;
          return tmp_;
        }
        return *stdstring_itr_;
      }
      pointer operator->() { return &(operator*()); }
      bool operator==(const self_type& rhs)
      {
        if (cstring_itr_)
          return cstring_itr_ == rhs.cstring_itr_;
        return stdstring_itr_ == rhs.stdstring_itr_;
      }
      bool operator!=(const self_type& rhs) { return !(*this == rhs); }
    private:
      const char** cstring_itr_;
      std::vector<std::string>::const_iterator stdstring_itr_;
      std::string tmp_;
    };

    virtual ~reader_base() {}

    operator bool() const
    {
      return this->good();
    }

    bool good() const
    {
      if (sav_reader_)
        return sav_reader_->good();
      else if (vcf_reader_)
        return vcf_reader_->good();
      return false;
    }

    bool fail() const
    {
      if (sav_reader_)
        return sav_reader_->fail();
      else if (vcf_reader_)
        return vcf_reader_->fail();
      return true;
    }

    bool bad() const
    {
      if (sav_reader_)
        return sav_reader_->bad();
      else if (vcf_reader_)
        return vcf_reader_->bad();
      return true;
    }

    template <typename T>
    bool read_variant(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN())
    {
      if (sav_reader_)
        sav_reader_->read_variant(destination, missing_value);
      else if (vcf_reader_)
        vcf_reader_->read_variant(destination, missing_value);

      return good();
    }

    std::vector<std::string> prop_fields() const;
    sample_iterator samples_begin() const;
    sample_iterator samples_end() const;
    std::size_t sample_size() const;
  protected:
    std::unique_ptr<sav::reader_base> sav_reader_;
    std::unique_ptr<vcf::reader_base> vcf_reader_;
  };

  class reader : public reader_base
  {
  public:
    reader() {}
    reader(const std::string& file_path);
    ~reader() {}

    template <typename T>
    reader& operator>>(T& destination);

    template <typename T>
    reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());
    std::vector<std::string> chromosomes() const;
  };

  class indexed_reader : public reader_base
  {
  public:
    indexed_reader() {}
    indexed_reader(const std::string& file_path, const region& reg);
    void reset_region(const region& reg);

    template <typename T>
    indexed_reader& operator>>(T& destination);

    template <typename T>
    indexed_reader& read(T& destination, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());

    template <typename T, typename Pred>
    indexed_reader& read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value = std::numeric_limits<typename T::vector_type::value_type>::quiet_NaN());
  private:

  };

  template <typename T>
  reader& reader::operator>>(T& destination)
  {
    read_variant(destination);
    return *this;
  }

  template <typename T>
  reader& reader::read(T& destination, const typename T::vector_type::value_type missing_value)
  {
    read_variant(destination, missing_value);
    return *this;
  }

  template <typename T>
  indexed_reader& indexed_reader::operator>>(T& destination)
  {
    read_variant(destination);
    return *this;
  }

  template <typename T>
  indexed_reader& indexed_reader::read(T& destination, const typename T::vector_type::value_type missing_value)
  {
    read_variant(destination, missing_value);
    return *this;
  }

  template <typename T, typename Pred>
  indexed_reader& indexed_reader::read_if(T& destination, Pred fn, const typename T::vector_type::value_type missing_value)
  {
    if (sav_reader_)
      dynamic_cast<sav::indexed_reader*>(sav_reader_.get())->read_if(destination, fn, missing_value);
    else if (vcf_reader_)
      dynamic_cast<vcf::indexed_reader*>(vcf_reader_.get())->read_if(destination, fn, missing_value);

    return *this;
  }

  template <typename VecType>
  using allele_variant_iterator =  basic_allele_variant_iterator<reader_base, VecType>;

  template <typename ValType>
  using dense_allele_variant_iterator =  basic_allele_variant_iterator<reader_base, std::vector<ValType>>;
  template <typename ValType>
  using sparse_allele_variant_iterator =  basic_allele_variant_iterator<reader_base, compressed_vector<ValType>>;
}

#endif //VC_READER_HPP