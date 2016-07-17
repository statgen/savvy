#ifndef VC_READER_HPP
#define VC_READER_HPP

#include <iterator>
#include <string>

namespace vc
{
  class variant
  {
  public:
  };

  class marker
  {
  public:
    virtual double calculate_allele_frequency() const = 0;
  };



  class reader
  {
  public:
    class const_iterator
    {
    public:
      typedef std::ptrdiff_t difference_type; //almost always ptrdif_t
      typedef marker value_type; //almost always T
      typedef const value_type& reference; //almost always T& or const T&
      typedef const value_type* pointer; //almost always T* or const T*
      typedef std::input_iterator_tag iterator_category;  //usually std::forward_iterator_tag or similar
    };

    reader(std::istream& input_stream);
    virtual ~reader() {}
  private:
  };



}
#endif //VC_READER_HPP
