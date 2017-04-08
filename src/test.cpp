
#include "cmf_reader.hpp"
#include "m3vcf_reader.hpp"
#include "vcf_reader.hpp"
#include "test_class.hpp"
#include "varint.hpp"
#include "vc.hpp"
#include "variant_iterator.hpp"
#include "reader.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <synced_bcf_reader.h>
#include <tbx.h>
#include <hts.h>

bool has_extension(const std::string& fullString, const std::string& ext)
{
  if (fullString.length() >= ext.length())
    return (0 == fullString.compare (fullString.length() - ext.length(), ext.length(), ext));
  else
    return false;
}

template <typename T>
void handle_file_reader(T& reader)
{
  typename T::input_iterator::buffer buff;
  typename T::input_iterator end;
  typename T::input_iterator it(reader, buff);

  while (it != end)
  {
    double af = it->calculate_allele_frequency();
    for (auto jt = it->begin(); jt != it->end(); ++jt)
    {
      vc::allele_status foo = *jt;
    }

    std::for_each(it->begin(), it->end(), [](const typename vc::allele_status& s)
    {
      vc::allele_status foo = s;
    });

    ++it;
  }

}

//int reader_tests()
//{
//  //----------------------------------------------------------------//
//  {
//    std::string file_path = "/foobar.cmf";
//    if (has_extension(file_path, ".cmf"))
//    {
//      vc::sav::reader input("/foobar.cmf");
//      handle_file_reader(input);
//    }
//    else if (has_extension(file_path, ".m3vcf"))
//    {
//      std::ifstream ifs("/foobar.m3vcf");
//      vc::m3vcf::reader input(ifs);
//      handle_file_reader(input);
//    }
//    else if (has_extension(file_path, ".vcf") || has_extension(file_path, "vcf.gz") || has_extension(file_path, ".bcf"))
//    {
//      vc::vcf::block buff;
//      vc::vcf::reader input(file_path);
//      vc::vcf::reader::input_iterator eof;
//      vc::vcf::reader::input_iterator cur(input, buff);
//
//      while (cur != eof)
//        ++cur;
//    }
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  {
//    std::string file_path = "/foobar.cmf";
//    if (has_extension(file_path, ".cmf"))
//    {
//      vc::sav::reader input("/foobar.cmf");
//      auto analysis = make_analysis(input);
//      analysis.run();
//    }
//    else
//    {
//      std::ifstream ifs("/foobar.m3vcf");
//      vc::m3vcf::reader input(ifs);
//      some_analysis<vc::m3vcf::reader> analysis(input);
//      analysis.run();
//    }
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  {
//    vc::sav::reader input("/foobar.cmf");
//    vc::sav::marker buff;
//
//    for (vc::sav::reader::input_iterator i(input, buff), eof; i != eof; ++i)
//    {
//      for (auto j = i->begin(); j != i->end(); ++j)
//      {
//
//      }
//    }
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  {
//    std::ifstream ifs("/foobar.m3vcf");
//    vc::m3vcf::reader input(ifs);
//    vc::m3vcf::block buff;
//
//    std::for_each(vc::m3vcf::reader::input_iterator(input, buff), vc::m3vcf::reader::input_iterator(), [](const vc::m3vcf::marker& m)
//    {
//      std::for_each(m.begin(), m.end(), [](const vc::allele_status& s)
//      {
//
//      });
//    });
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  vc::sav::marker m;
//  std::uint64_t ploidy_level = 2;
//  std::uint64_t sample_size = 1000;
//  std::vector<int> zero_one_two_vec(sample_size, 0);
//
//  std::for_each(m.non_ref_begin(), m.non_ref_end(), [&zero_one_two_vec, ploidy_level](const vc::sav::marker::sparse_vector_allele& a)
//  {
//    if (a.status == vc::allele_status::has_alt)
//      ++(zero_one_two_vec[a.offset / ploidy_level]);
//  });
//
//  std::size_t i = 0;
//  std::for_each(m.begin(), m.end(), [&zero_one_two_vec, &i, ploidy_level](const vc::allele_status& s)
//  {
//    if (s == vc::allele_status::has_alt)
//      ++(zero_one_two_vec[i / ploidy_level]);
//  });
//
//  i = 0;
//  for (const auto& s : m)
//  {
//    if (s == vc::allele_status::has_alt)
//      ++(zero_one_two_vec[i / ploidy_level]);
//  }
//  //----------------------------------------------------------------//
//  return 0;
//}

int varint_test()
{
  std::vector<std::uint64_t> arr(0xFFFFFF);
  for (std::uint64_t i = 0; i < arr.size(); ++i)
    arr[i] = i;
  std::cout << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;

  {
    std::ofstream non_compressed_arr_ostream("foo-not-compressed.bin", std::ios::binary);
    const auto encode_start = std::chrono::high_resolution_clock::now();
    non_compressed_arr_ostream.write((char*)arr.data(), arr.size() * sizeof(std::uint64_t));
    non_compressed_arr_ostream.flush();
    //std::copy(arr.begin(), arr.end(), back_it);
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    non_compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream non_compressed_arr_istream("foo-not-compressed.bin", std::ios::binary);
    const auto decode_start = std::chrono::high_resolution_clock::now();
    non_compressed_arr_istream.read((char*)arr.data(), arr.size() * sizeof(std::uint64_t));
//    for (std::size_t i = 0; i < arr.size(); ++i)
//      arr[i] = ntohll(arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "Non-compressed copy: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    std::cout << std::endl;
  }

  {
    std::ofstream compressed_arr_ostream("foo-0bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::varint_encode(i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-0bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::varint_decode(decode_it, std::istreambuf_iterator<char>(), arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "0-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-1bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::one_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-1bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::one_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "1-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-2bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::two_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-2bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::two_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "2-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-3bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::three_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-3bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::three_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "3-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-4bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::four_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-4bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::four_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "4-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-5bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::five_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-5bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::five_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "5-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-6bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::six_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-6bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::six_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "6-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  {
    std::uint8_t prefix_data = 0;
    std::ofstream compressed_arr_ostream("foo-7bit.bin");
    const auto encode_start = std::chrono::high_resolution_clock::now();
    std::ostreambuf_iterator<char> output_it(compressed_arr_ostream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      vc::seven_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-7bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++vc::seven_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
    auto decode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - decode_start).count();
    std::cout << "7-bit prefixed: " << std::accumulate(arr.begin(), arr.end(), 0ULL) << std::endl;
    std::cout << "Decode elapsed time: " << decode_elapsed_time << "ms" << std::endl;
    compressed_arr_istream.close();
    std::cout << std::endl;
  }

  return 0;
}

template <typename Proc>
class timed_procedure_call
{
public:
  timed_procedure_call(Proc& procedure)
  {
    start_ = std::chrono::high_resolution_clock::now();
    return_value_ = procedure();
    end_ = std::chrono::high_resolution_clock::now();
  }
  bool return_value() const { return return_value_; }
  template <typename Duration>
  long long int elapsed_time()
  {
    return std::chrono::duration_cast<Duration>(end_ - start_).count();
  }
private:
  bool return_value_;
  std::chrono::high_resolution_clock::time_point start_;
  std::chrono::high_resolution_clock::time_point end_;
};

template <typename Proc>
timed_procedure_call<Proc> time_procedure(Proc& p)
{
  return timed_procedure_call<Proc>(p);
}

template <typename R1, typename R2>
class file_checksum_test
{
public:
  file_checksum_test(R1& reader1, R2& reader2) : reader1_(reader1), reader2_(reader2) {}
  bool operator()() const
  {
    std::size_t checksum1 = get_checksum(reader1_);
    std::size_t checksum2 = get_checksum(reader2_);

    std::cout << checksum1 << " " << checksum2 << std::endl;

    return checksum1 == checksum2;
  }
private:
  template <typename T>
  static std::size_t hash_combine(std::size_t seed, const T& val)
  {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;
  }

  template <typename ReaderType>
  static std::size_t get_checksum(ReaderType& reader)
  {
    std::size_t ret = 0;

    vc::basic_variant_iterator<ReaderType, std::vector<float>> cur(reader);
    vc::basic_variant_iterator<ReaderType, std::vector<float>> end;

    std::size_t num_markers = 0;
    while (cur != end)
    {
      ret = hash_combine(ret, cur->locus());
      ret = hash_combine(ret, cur->ref());
      ret = hash_combine(ret, cur->alt());

      for (auto gt = cur->begin(); gt != cur->end(); ++gt)
        ret = hash_combine(ret, static_cast<int>(*gt));

      ++cur;
      ++num_markers;
    }
    std::cout << "Marker Count: " << num_markers << std::endl;

    return ret;
  }
  R1& reader1_;
  R2& reader2_;
};

template <typename T1, typename T2>
file_checksum_test<T1, T2> make_file_checksum_test(T1& a, T2& b)
{
  return file_checksum_test<T1, T2>(a, b);
}

void run_file_checksum_test()
{

  vc::open_files(std::make_tuple("test_file.vcf", "test_file.sav"), [](auto&& input_file_reader1, auto&& input_file_reader2)
  {
    auto t = make_file_checksum_test(input_file_reader1, input_file_reader2);
    std::cout << "Starting checksum test ..." << std::endl;
    auto timed_call = time_procedure(t);
    std::cout << "Returned: " << (timed_call.return_value() ? "True" : "FALSE") << std::endl;
    std::cout << "Elapsed Time: " << timed_call.template elapsed_time<std::chrono::milliseconds>() << "ms" << std::endl;
  });
}

void convert_file_test()
{
  {
    vc::vcf::reader input("test_file.vcf");
    vc::vcf::variant_iterator<std::vector<float>> cur(input);
    vc::vcf::variant_iterator<std::vector<float>> eof;
    const std::string chrom = cur != eof ? cur->chromosome() : "";
    const std::size_t ploidy = cur != eof ? (cur->size() / input.sample_count()) : 0;

    vc::sav::writer compact_output("test_file.sav", chrom, ploidy, input.samples_begin(), input.samples_end(), input.prop_fields_begin(), input.prop_fields_end());

    while (cur != eof)
    {
      if (cur->chromosome() != chrom)
      {
        std::cerr << "Multiple chromosomes encountered. CMF files can only contain one chromosome." << std::endl;
        break;
      }

      //    vc::haplotype_vector<std::vector<float>> m(std::string(chrom), cur->locus(), std::string(cur->ref()), std::string(cur->alt()), sample_ids.size(), ploidy, std::vector<float>());
      //    for (auto it = cur->begin(); it != cur->end(); ++it)
      //    {
      //      switch (*it)
      //      {
      //        case vc::allele_status::has_alt:
      //          m[std::distance(cur->begin(), it)] = 1.0;
      //          break;
      //        case vc::allele_status::is_missing:
      //          m[std::distance(cur->begin(), it)] = std::numeric_limits<float>::quiet_NaN();
      //          break;
      //      }
      //    }

      compact_output << *cur;

      ++cur;
    }
  }

  run_file_checksum_test();
}

//template <typename M>
//double inner_product(const M& mrkr, std::vector<double>& vec, const double start_val = 0.0)
//{
//  double ret = start_val;
//
//  std::size_t i = 0;
//  for (const vc::allele_status& gt : mrkr)
//  {
//    if (gt == vc::allele_status::has_alt)
//      ret += 1.0 * vec[i];
//    ++i;
//  }
//
//  return ret;
//}
//
//template <>
//double inner_product<vc::sav::marker>(const vc::sav::marker& mrkr, std::vector<double>& vec, const double start_val)
//{
//  double ret = start_val;
//
//  for (auto it = mrkr.non_ref_begin(); it != mrkr.non_ref_end(); ++it)
//  {
//    if (it->status == vc::allele_status::has_alt)
//      ret += 1.0 * vec[it->offset];
//    else
//      ret += 0.04 * vec[it->offset];
//  }
//
//  return ret;
//}

//class file_handler_functor
//{
//public:
//  template <typename T>
//  void operator()(T&& input_file_reader)
//  {
//    typedef typename T::input_iterator input_iterator_type;
//    typename T::input_iterator::buffer buf{};
//
//    for (auto it = input_iterator_type(input_file_reader, buf); it != input_iterator_type(); ++it)
//    {
//      inner_product(*it, phenotypes_);
//    }
//  };
//private:
//  std::vector<double> phenotypes_;
//};

class triple_file_handler_functor
{
public:
  template <typename T, typename T2, typename T3>
  void operator()(T&& input_file_reader, T2&& input_file_reader2, T3&& input_file_reader3)
  {
    input_file_reader.sample_count();
    input_file_reader2.sample_count();
    input_file_reader3.sample_count();

    std::tuple<T, T2, T3> file_readers(std::move(input_file_reader), std::move(input_file_reader2), std::move(input_file_reader3));

  }
};

class marker_handler_functor
{
public:
  template <typename T>
  void operator()(const T& mrkr)
  {
    std::uint64_t pos = mrkr.pos();
    std::string ref = mrkr.ref();
    std::string alt = mrkr.alt();
    std::for_each(mrkr.begin(), mrkr.end(), [](const vc::allele_status& s)
    {

    });
  }
};

//class marker_counter
//{
//public:
//  marker_counter() = default;
//  marker_counter(const marker_counter&) = delete;
//  marker_counter(marker_counter&&) = delete;
//  marker_counter& operator=(const marker_counter&) = delete;
//  marker_counter& operator=(marker_counter&&) = delete;
//  template <typename T, typename T2>
//  void operator()(T&& input_file_reader, T2&& input_file_reader2)
//  {
//    typename T::input_iterator::buffer buff;
//    typename T::input_iterator cur(input_file_reader, buff);
//    typename T::input_iterator end;
//
//    typename T2::input_iterator::buffer buff2;
//    typename T2::input_iterator cur2(input_file_reader2, buff2);
//    typename T2::input_iterator end2;
//
//    inner_product(*cur2, *end2);
//    while (cur != end)
//    {
//      ++file1_cnt_;
//      ++cur;
//    }
//
//    while (cur2 != end2)
//    {
//      ++file2_cnt_;
//      ++cur2;
//    }
//
//  }
//  std::size_t file1_count() const { return file1_cnt_; }
//  std::size_t file2_count() const { return file2_cnt_; }
//private:
//  std::size_t file1_cnt_ = 0;
//  std::size_t file2_cnt_ = 0;
//};

void random_access_test()
{
  vc::sav::writer::create_index("test_file.sav");

  vc::indexed_reader rdr("test_file.sav", {"20", 17000, 1120000});
  auto it = vc::variant_iterator<std::vector<float>>(rdr);
  auto end = vc::variant_iterator<std::vector<float>>{};
  for ( ; it != end; ++it)
  {
    std::cout << it->locus() << std::endl;
  }

  rdr.reset_region({"20", 17000, 1120000});
  vc::dense_haplotype_vector<float> v;
  while (rdr >> v)
  {
    std::cout << v.locus() << std::endl;
  }
}

void generic_reader_test()
{
  vc::reader rdr1("test_file.sav");
  vc::reader rdr2("test_file.vcf");


  auto t = make_file_checksum_test(rdr1, rdr2);
  std::cout << "Starting checksum test ..." << std::endl;
  auto timed_call = time_procedure(t);
  std::cout << "Returned: " << (timed_call.return_value() ? "True" : "FALSE") << std::endl;
  std::cout << "Elapsed Time: " << timed_call.template elapsed_time<std::chrono::milliseconds>() << "ms" << std::endl;
}

//#include <boost/numeric/ublas/vector_sparse.hpp>

int main(int argc, char** argv)
{
//  using namespace boost::numeric;
//  ublas::compressed_vector<float> v(10);
//  v[3] = 0;
//  v[5] = 3;
//  std::cout << v.nnz() << std::endl;
//  vc::haplotype_vector<ublas::compressed_vector<float>> v2;
//
//  std::cout << std::endl;
//  return 0;
  std::cout << "[0] Run all tests." << std::endl;
  std::cout << "[1] Run varint test." << std::endl;
  std::cout << "[2] Run file conversion test." << std::endl;
  std::cout << "[3] Run generic reader test." << std::endl;
  std::cout << "[4] Run random access test." << std::endl;

  char t = '0';
  std::cin >> t;

  switch (t)
  {
    case '0':
      varint_test();
      convert_file_test();
      break;
    case '1':
      varint_test();
      break;
    case '2':
      convert_file_test();
      break;
    case '3':
      generic_reader_test();
      break;
    case '4':
      random_access_test();
      break;
    default:
      std::cout << "Invalid Input" << std::endl;
  }



//  vc::open_marker_files(triple_file_handler_functor(), "chr1.bcf", "chr1.cmf", "chr1.m3vcf");
//
//  vc::open_marker_files(std::make_tuple("chr1.cmf", "chr1.m3vcf"), [](auto&& input_file_reader1, auto&& input_file_reader2)
//  {
//    typedef typename std::remove_reference<decltype(input_file_reader1)>::type R1;
//    typename R1::input_iterator::buffer buf{};
//    typename R1::input_iterator eof{};
//    typename R1::input_iterator it(input_file_reader1, buf);
//
//    typedef typename std::remove_reference<decltype(input_file_reader2)>::type R2;
//    typename R2::input_iterator::buffer buf2{};
//    typename R2::input_iterator eof2{};
//    typename R2::input_iterator it2(input_file_reader2, buf2);
//
//    while (it != eof)
//    {
//
//      ++it;
//    }
//
//    while (it2 != eof2)
//    {
//
//      ++it2;
//    }
//
//  });
//
//  vc::open_marker_file("chr1.bcf", [](auto&& input_file_reader)
//  {
//    typedef typename std::remove_reference<decltype(input_file_reader)>::type R;
//    typename R::input_iterator::buffer buf{};
//    typename R::input_iterator eof{};
//    typename R::input_iterator it(input_file_reader, buf);
//
//    while (it != eof)
//    {
//      it->pos();
//      it->ref() + ":" + it->alt();
//      for (const vc::allele_status& gt : *it)
//      {
//
//      }
//      ++it;
//    }
//  });
//
//  vc::open_marker_file("chr1.bcf", file_handler_functor());
//  file_handler_functor f;
//  vc::open_marker_file("chr1.bcf", f);
//
//  vc::iterate_marker_file("chr1.bcf", marker_handler_functor());



  return 0;
}

