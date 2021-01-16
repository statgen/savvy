/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/sav_reader.hpp"
#include "savvy/m3vcf_reader.hpp"
#include "savvy/vcf_reader.hpp"
#include "test/test_class.hpp"
#include "savvy/varint.hpp"
#include "savvy/savvy.hpp"
#include "savvy/variant_iterator.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"
#include "savvy/site_info.hpp"
#include "savvy/data_format.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <sys/stat.h>


//bool has_extension(const std::string& fullString, const std::string& ext)
//{
//  if (fullString.length() >= ext.length())
//    return (0 == fullString.compare (fullString.length() - ext.length(), ext.length(), ext));
//  else
//    return false;
//}
//
//template <typename T>
//void handle_file_reader(T& reader)
//{
//  typename T::input_iterator::buffer buff;
//  typename T::input_iterator end;
//  typename T::input_iterator it(reader, buff);
//
//  while (it != end)
//  {
//    double af = it->calculate_allele_frequency();
//    for (auto jt = it->begin(); jt != it->end(); ++jt)
//    {
//      savvy::allele_status foo = *jt;
//    }
//
//    std::for_each(it->begin(), it->end(), [](const typename savvy::allele_status& s)
//    {
//      savvy::allele_status foo = s;
//    });
//
//    ++it;
//  }
//
//}

//int reader_tests()
//{
//  //----------------------------------------------------------------//
//  {
//    std::string file_path = "/foobar.cmf";
//    if (has_extension(file_path, ".cmf"))
//    {
//      savvy::sav::reader input("/foobar.cmf");
//      handle_file_reader(input);
//    }
//    else if (has_extension(file_path, ".m3vcf"))
//    {
//      std::ifstream ifs("/foobar.m3vcf");
//      savvy::m3vcf::reader input(ifs);
//      handle_file_reader(input);
//    }
//    else if (has_extension(file_path, ".vcf") || has_extension(file_path, "vcf.gz") || has_extension(file_path, ".bcf"))
//    {
//      savvy::vcf::block buff;
//      savvy::vcf::reader input(file_path);
//      savvy::vcf::reader::input_iterator eof;
//      savvy::vcf::reader::input_iterator cur(input, buff);
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
//      savvy::sav::reader input("/foobar.cmf");
//      auto analysis = make_analysis(input);
//      analysis.run();
//    }
//    else
//    {
//      std::ifstream ifs("/foobar.m3vcf");
//      savvy::m3vcf::reader input(ifs);
//      some_analysis<savvy::m3vcf::reader> analysis(input);
//      analysis.run();
//    }
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  {
//    savvy::sav::reader input("/foobar.cmf");
//    savvy::sav::marker buff;
//
//    for (savvy::sav::reader::input_iterator i(input, buff), eof; i != eof; ++i)
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
//    savvy::m3vcf::reader input(ifs);
//    savvy::m3vcf::block buff;
//
//    std::for_each(savvy::m3vcf::reader::input_iterator(input, buff), savvy::m3vcf::reader::input_iterator(), [](const savvy::m3vcf::marker& m)
//    {
//      std::for_each(m.begin(), m.end(), [](const savvy::allele_status& s)
//      {
//
//      });
//    });
//  }
//  //----------------------------------------------------------------//
//
//  //----------------------------------------------------------------//
//  savvy::sav::marker m;
//  std::uint64_t ploidy_level = 2;
//  std::uint64_t sample_size = 1000;
//  std::vector<int> zero_one_two_vec(sample_size, 0);
//
//  std::for_each(m.non_ref_begin(), m.non_ref_end(), [&zero_one_two_vec, ploidy_level](const savvy::sav::marker::sparse_vector_allele& a)
//  {
//    if (a.status == savvy::allele_status::has_alt)
//      ++(zero_one_two_vec[a.offset / ploidy_level]);
//  });
//
//  std::size_t i = 0;
//  std::for_each(m.begin(), m.end(), [&zero_one_two_vec, &i, ploidy_level](const savvy::allele_status& s)
//  {
//    if (s == savvy::allele_status::has_alt)
//      ++(zero_one_two_vec[i / ploidy_level]);
//  });
//
//  i = 0;
//  for (const auto& s : m)
//  {
//    if (s == savvy::allele_status::has_alt)
//      ++(zero_one_two_vec[i / ploidy_level]);
//  }
//  //----------------------------------------------------------------//
//  return 0;
//}

bool file_exists(const std::string& file_path)
{
  struct stat st;
  return (stat(file_path.c_str(), &st) == 0);
}

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
      savvy::varint_encode(i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-0bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::varint_decode(decode_it, std::istreambuf_iterator<char>(), arr[i]);
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
      savvy::one_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-1bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::one_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::two_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-2bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::two_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::three_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-3bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::three_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::four_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-4bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::four_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::five_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-5bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::five_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::six_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-6bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::six_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
      savvy::seven_bit_prefixed_varint::encode(prefix_data, i, output_it);
    compressed_arr_ostream.flush();
    auto encode_elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - encode_start).count();
    std::cout << "Encode elapsed time: " << encode_elapsed_time << "ms" << std::endl;
    compressed_arr_ostream.close();

    std::fill(arr.begin(), arr.end(), 0);

    std::ifstream compressed_arr_istream("foo-7bit.bin");
    const auto decode_start = std::chrono::high_resolution_clock::now();
    std::istreambuf_iterator<char> decode_it(compressed_arr_istream);
    for (std::uint64_t i = 0; i < arr.size(); ++i)
      decode_it = ++savvy::seven_bit_prefixed_varint::decode(decode_it, std::istreambuf_iterator<char>(), prefix_data, arr[i]);
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
  file_checksum_test(R1& reader1, R2& reader2, const std::string& fmt) : reader1_(reader1), reader2_(reader2), fmt_field_(fmt) {}
  bool operator()() const
  {
    std::size_t checksum1 = get_checksum(reader1_, fmt_field_);
    std::size_t checksum2 = get_checksum(reader2_, fmt_field_);

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
  static std::size_t get_checksum(ReaderType& reader, const std::string& fmt_field)
  {
    std::size_t ret = 0;

    savvy::variant var;
    std::vector<float> data;

    //auto prop_fields = reader.info_fields();

    //prop_fields.erase(std::find(prop_fields.begin(), prop_fields.end(), "AF"));

    std::size_t num_markers = 0;
    while (reader.read(var))
    {
      ret = hash_combine(ret, var.position());
      ret = hash_combine(ret, var.ref());
      for (const auto& a : var.alts())
        ret = hash_combine(ret, a);

      for (const auto& field : var.info())
      {
        std::stringstream prop_val;
        prop_val << field.second;
        std::string prop_val_str = prop_val.str();
        if (field.first == "AF")
          prop_val_str = prop_val_str.substr(0, std::min((std::size_t)3, prop_val_str.size()));

        ret = hash_combine(ret, prop_val_str);
      }

      var.get_format(fmt_field, data);

      for (auto gt = data.begin(); gt != data.end(); ++gt)
        ret = hash_combine(ret, savvy::sav::detail::allele_encoder<7>::encode(*gt));

      //print_variant(prop_fields, anno, data);

      ++num_markers;
    }
    std::cout << "Marker Count: " << num_markers << std::endl;

    return ret;
  }
  R1& reader1_;
  R2& reader2_;
  std::string fmt_field_;
};

template <typename T1, typename T2>
file_checksum_test<T1, T2> make_file_checksum_test(T1& a, T2& b, const std::string& fmt)
{
  return file_checksum_test<T1, T2>(a, b, fmt);
}

void run_file_checksum_test(const std::string f1, const std::string f2, const std::string& fmt)
{
  savvy::reader input_file_reader1(f1);
  //input_file_reader1.set_policy(savvy::vcf::empty_vector_policy::skip);
  savvy::reader input_file_reader2(f2);
  auto t = make_file_checksum_test(input_file_reader1, input_file_reader2, fmt);
  std::cout << "Starting checksum test ..." << std::endl;
  auto timed_call = time_procedure(t);
  std::cout << "Returned: " << (timed_call.return_value() ? "True" : "FALSE") << std::endl;
  std::cout << "Elapsed Time: " << timed_call.template elapsed_time<std::chrono::milliseconds>() << "ms" << std::endl;
  assert(timed_call.return_value());
}


void convert_file_test(const std::string& fmt_field)
{
  {
    savvy::reader input(SAVVYT_VCF_FILE);

    savvy::variant var;
    savvy::compressed_vector<float> data;

    auto file_info = input.headers();
    file_info.reserve(file_info.size() + 3);
    file_info.insert(file_info.begin(), {"INFO", "<ID=FILTER,Description=\"Variant filter\">"});
    file_info.insert(file_info.begin(), {"INFO", "<ID=QUAL,Description=\"Variant quality\">"});
    file_info.insert(file_info.begin(), {"INFO", "<ID=ID,Description=\"Variant ID\">"});
    savvy::writer output(fmt_field == "HDS" ? SAVVYT_SAV_FILE_DOSE : SAVVYT_SAV_FILE_HARD, savvy::file::format::sav2, file_info, input.samples());

    std::size_t cnt = 0;
    while (input.read(var))
    {
      output.write(var);
      ++cnt;
    }

    assert(output.good() && !input.bad());
    assert(cnt == SAVVYT_MARKER_COUNT_HARD);
  }

  run_file_checksum_test(SAVVYT_VCF_FILE, fmt_field == "HDS" ? SAVVYT_SAV_FILE_DOSE : SAVVYT_SAV_FILE_HARD, fmt_field);

  //savvy::sav::writer::create_index(fmt_field == "HDS" ? SAVVYT_SAV_FILE_DOSE : SAVVYT_SAV_FILE_HARD);
}

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


void sav_random_access_test(const std::string& fmt_field)
{
  savvy::reader rdr(fmt_field == "GT" ? SAVVYT_SAV_FILE_HARD : SAVVYT_SAV_FILE_DOSE);
  assert(rdr.good());
  rdr.reset_bounds({"20", 1234600, 2234567});
  assert(rdr.good());

  savvy::variant anno;
  std::vector<float> buf;

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "20");
  assert(anno.position() == 1234667);
  assert(anno.ref() == "G");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "A");

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "20");
  assert(anno.position() == 1234767);
  assert(anno.ref() == "T");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "A");

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "20");
  assert(anno.position() == 2230237);
  assert(anno.ref() == "T");
  assert(anno.alts().size() == 0);


  assert(rdr.read(anno));
  assert(fmt_field != "GT" || anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "20");
  assert(anno.position() == 2234567);
  assert(anno.ref() == "GTC");
  assert(anno.alts().size() == 2);
  assert(anno.alts()[0] == "G");
  assert(anno.alts()[1] == "GTCT");

  assert(rdr.good());

  assert(!rdr.read(anno));

  rdr.reset_bounds({"18", 2234600, 2234700});

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "18");
  assert(anno.position() == 2234668);
  assert(anno.ref() == "G");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "A");

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "18");
  assert(anno.position() == 2234679);
  assert(anno.ref() == "G");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "T");

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "18");
  assert(anno.position() == 2234687);
  assert(anno.ref() == "G");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "A");

  assert(rdr.read(anno));
  assert(anno.get_format(fmt_field, buf));
  assert(anno.chromosome() == "18");
  assert(anno.position() == 2234697);
  assert(anno.ref() == "T");
  assert(anno.alts().size() == 1);
  assert(anno.alts()[0] == "A");

  assert(rdr.good());

  assert(!rdr.read(anno));
}

void generic_reader_test(const std::string& path, const std::string& fmt_field, std::size_t expected_markers)
{
  savvy::reader rdr(path);
  assert(rdr.good());

  std::vector<std::string> subset = {"NA00003","NA00005", "FAKE_ID"};
  auto intersect = rdr.subset_samples({subset.begin(), subset.end()});
  assert(intersect.size() == 2);

  std::size_t stride = 2;
  savvy::variant i;
  std::vector<float> d;
  std::size_t cnt{};
  while (rdr.read(i))
  {
    i.get_format(fmt_field, d);
    assert(d.size() == intersect.size() * stride);
    ++cnt;
  }
  assert(cnt == expected_markers);
}

template <typename R>
void subset_test(const std::string& path, const std::string& fmt_field)
{
  R rdr(path);
  assert(rdr.good());

  std::vector<std::string> subset = {"NA00003","NA00005", "FAKE_ID"};
  auto intersect = rdr.subset_samples({subset.begin(), subset.end()});
  assert(intersect.size() == 2);

  savvy::variant i;
  std::vector<float> d;
  std::size_t cnt{};
  while (rdr.read(i))
  {
    i.get_format(fmt_field, d);
    assert(d.size() == intersect.size() * 2);
    ++cnt;
  }
  assert(cnt == SAVVYT_MARKER_COUNT_HARD);
}

int main(int argc, char** argv)
{
  std::string cmd = (argc < 2) ? "" : argv[1];

  if (cmd.empty())
  {
    std::cout << "Enter Command:" << std::endl;
    std::cout << "- convert-file" << std::endl;
    std::cout << "- generic-reader" << std::endl;
    std::cout << "- random-access" << std::endl;
    std::cout << "- subset" << std::endl;
    std::cout << "- varint" << std::endl;
    std::cin >> cmd;
  }


  if (cmd == "convert-file")
  {
    convert_file_test("GT");
    convert_file_test("HDS");
  }
//  else if (cmd == "generic-reader")
//  {
//    if (!file_exists(SAVVYT_SAV_FILE_HARD)) convert_file_test("GT");
//    if (!file_exists(SAVVYT_SAV_FILE_DOSE)) convert_file_test("HDS");
//
//    generic_reader_test(SAVVYT_VCF_FILE, "GT", SAVVYT_MARKER_COUNT_HARD);
//    generic_reader_test(SAVVYT_VCF_FILE, "HDS", SAVVYT_MARKER_COUNT_DOSE);
//
//    generic_reader_test(SAVVYT_SAV_FILE_HARD, "GT", SAVVYT_MARKER_COUNT_HARD);
//    //generic_reader_test(SAVVYT_SAV_FILE_HARD, "HDS", SAVVYT_MARKER_COUNT_HARD);
//
//    //generic_reader_test(SAVVYT_SAV_FILE_DOSE, "GT", SAVVYT_MARKER_COUNT_DOSE);
//    generic_reader_test(SAVVYT_SAV_FILE_DOSE, "HDS", SAVVYT_MARKER_COUNT_DOSE);
//  }
  else if (cmd == "random-access")
  {
    if (!file_exists(SAVVYT_SAV_FILE_HARD)) convert_file_test("GT");
    if (!file_exists(SAVVYT_SAV_FILE_DOSE)) convert_file_test("HDS");

    sav_random_access_test("GT");
    sav_random_access_test("HDS");
  }
  else if (cmd == "subset")
  {
    if (!file_exists(SAVVYT_SAV_FILE_HARD)) convert_file_test("GT");

    subset_test<savvy::reader>(SAVVYT_VCF_FILE, "GT");
    subset_test<savvy::reader>(SAVVYT_SAV_FILE_HARD, "GT");
  }
  else if (cmd == "varint")
  {
    varint_test();
  }
  else
  {
    std::cerr << "Invalid Command" << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}

