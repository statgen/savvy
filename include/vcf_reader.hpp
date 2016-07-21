#ifndef LIBVC_VCF_READER_HPP
#define LIBVC_VCF_READER_HPP

#include <iterator>
#include <string>


namespace vc
{
namespace vcf
{
#include "vcf.h"
}
}

namespace vc
{
  namespace vcf
  {

    class marker
    {
    public:
      marker();
      ~marker();

      static bool read_marker(marker& destination, htsFile* hts_file_, bcf_hdr_t* hts_hdr_);
    private:
      bcf1_t* hts_rec_;
    };

    void foo();
    class reader
    {
    public:
      class input_iterator
      {
      public:
        typedef input_iterator self_type;
        typedef std::ptrdiff_t difference_type;
        typedef marker value_type;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::input_iterator_tag iterator_category;
        typedef marker buffer;

        input_iterator() : file_reader_(nullptr), buffer_(nullptr) {}
        input_iterator(reader& file_reader, marker& buffer) : file_reader_(&file_reader), buffer_(&buffer) {}
        void increment()
        {
          if (!file_reader_->read_next_marker(*buffer_))
            file_reader_ = nullptr;
        }
        self_type& operator++(){ increment(); return *this; }
        void operator++(int) { increment(); }
        reference operator*() { return *buffer_; }
        pointer operator->() { return buffer_; }
        bool operator==(const self_type& rhs) { return (file_reader_ == rhs.file_reader_); }
        bool operator!=(const self_type& rhs) { return (file_reader_ != rhs.file_reader_); }
      private:
        reader* file_reader_;
        marker* buffer_;
      };

      reader(const std::string& file_path);
      ~reader();
      bool read_next_marker(marker& destination);
    private:
      htsFile* hts_file_;
      bcf_hdr_t* hts_hdr_;
    };
  }
}

#endif //LIBVC_VCF_READER_HPP
