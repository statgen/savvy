
#ifndef LIBVC_VARIANT_ITERATOR_HPP
#define LIBVC_VARIANT_ITERATOR_HPP

#include "haplotype_vector.hpp"

#include <iterator>

namespace vc
{
  template <typename ReaderType, typename VectorType>
  class basic_variant_iterator
  {
  public:
    typedef basic_variant_iterator self_type;
    typedef std::ptrdiff_t difference_type;
    typedef haplotype_vector<VectorType> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef std::input_iterator_tag iterator_category;

    basic_variant_iterator() : file_reader_(nullptr) {}
    basic_variant_iterator(ReaderType& file_reader) :
      file_reader_(&file_reader)
    {
      increment();
    }

    void increment()
    {
      bool b = file_reader_->good();
      file_reader_->read(m_);
      if (!file_reader_->good())
        file_reader_ = nullptr;
    }
    self_type& operator++(){ increment(); return *this; }
    void operator++(int) { increment(); }
    reference operator*() { return m_; }
    pointer operator->() { return &m_; }
    bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
    bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
  private:
    ReaderType* file_reader_;
    value_type m_;
  };
}

#endif //LIBVC_VARIANT_ITERATOR_HPP
