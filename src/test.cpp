
#include "cvcf_reader.hpp"
#include "m3vcf_reader.hpp"
#include "vcf_reader.hpp"
#include "test_class.hpp"
#include "varint.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <sstream>

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

int reader_tests()
{
  //----------------------------------------------------------------//
  {
    std::string file_path = "/foobar.cvcf";
    if (has_extension(file_path, ".cvcf"))
    {
      std::ifstream ifs("/foobar.cvcf");
      vc::cvcf::reader input(ifs);
      handle_file_reader(input);
    }
    else if (has_extension(file_path, ".m3vcf"))
    {
      std::ifstream ifs("/foobar.m3vcf");
      vc::m3vcf::reader input(ifs);
      handle_file_reader(input);
    }
    else if (has_extension(file_path, ".vcf") || has_extension(file_path, "vcf.gz") || has_extension(file_path, ".bcf"))
    {
      vc::vcf::block buff;
      vc::vcf::reader input(file_path);
      vc::vcf::reader::input_iterator eof;
      vc::vcf::reader::input_iterator cur(input, buff);

      while (cur != eof)
        ++cur;
    }
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  {
    std::string file_path = "/foobar.cvcf";
    if (has_extension(file_path, ".cvcf"))
    {
      std::ifstream ifs("/foobar.cvcf");
      vc::cvcf::reader input(ifs);
      auto analysis = make_analysis(input);
      analysis.run();
    }
    else
    {
      std::ifstream ifs("/foobar.m3vcf");
      vc::m3vcf::reader input(ifs);
      some_analysis<vc::m3vcf::reader> analysis(input);
      analysis.run();
    }
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  {
    std::ifstream ifs("/foobar.cvcf");
    vc::cvcf::reader input(ifs);
    vc::cvcf::marker buff;

    for (vc::cvcf::reader::input_iterator i(input, buff), eof; i != eof; ++i)
    {
      for (auto j = i->begin(); j != i->end(); ++j)
      {

      }
    }
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  {
    std::ifstream ifs("/foobar.m3vcf");
    vc::m3vcf::reader input(ifs);
    vc::m3vcf::block buff;

    std::for_each(vc::m3vcf::reader::input_iterator(input, buff), vc::m3vcf::reader::input_iterator(), [](const vc::m3vcf::marker& m)
    {
      std::for_each(m.begin(), m.end(), [](const vc::allele_status& s)
      {

      });
    });
  }
  //----------------------------------------------------------------//

  //----------------------------------------------------------------//
  vc::cvcf::marker m;
  std::uint64_t ploidy_level = 2;
  std::uint64_t sample_size = 1000;
  std::vector<int> zero_one_two_vec(sample_size, 0);

  std::for_each(m.non_ref_begin(), m.non_ref_end(), [&zero_one_two_vec, ploidy_level](const vc::cvcf::marker::sparse_vector_allele& a)
  {
    if (a.status == vc::allele_status::has_alt)
      ++(zero_one_two_vec[a.offset / ploidy_level]);
  });

  std::size_t i = 0;
  std::for_each(m.begin(), m.end(), [&zero_one_two_vec, &i, ploidy_level](const vc::allele_status& s)
  {
    if (s == vc::allele_status::has_alt)
      ++(zero_one_two_vec[i / ploidy_level]);
  });

  i = 0;
  for (const auto& s : m)
  {
    if (s == vc::allele_status::has_alt)
      ++(zero_one_two_vec[i / ploidy_level]);
  }
  //----------------------------------------------------------------//
  return 0;
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

void convert_file_test()
{
  {
    std::ofstream ofs("sample_conversion.cvcf", std::ios::binary);

    vc::vcf::block buff;
    vc::vcf::reader input("sample_file.vcf");
    vc::vcf::reader::input_iterator eof;
    vc::vcf::reader::input_iterator cur(input, buff);

    std::cout << "CHROM\tPOS\tREF\tALT\t";
    for (auto it = input.samples_begin(); it != input.samples_end(); ++it)
    {
      std::cout << *it;
      if (it + 1 < input.samples_end())
        std::cout << "\t";
    }
    std::cout << std::endl;

    std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
    std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
    vc::cvcf::writer compact_output(ofs, "22", 2, sample_ids.begin(), sample_ids.end());

    while (cur != eof)
    {
      std::cout << vc::vcf::reader::get_chromosome(input, *cur) << "\t" << cur->pos() << "\t" << cur->ref() << "\t" << cur->alt() << "\t";
      for (auto gt = cur->begin(); gt != cur->end(); )
      {
        if (*gt == vc::allele_status::has_ref)
          std::cout << "0";
        else if (*gt == vc::allele_status::has_alt)
          std::cout << "1";
        else
          std::cout << ".";

        ++gt;
        if (gt != cur->end())
          std::cout << "\t";
      }

      std::cout << std::endl;

      compact_output << vc::cvcf::marker(cur->pos(), "", cur->ref(), cur->alt(), cur->begin(), cur->end());

      ++cur;
    }
  }

  {
    std::cout << std::endl << std::endl;
    std::ifstream ifs("sample_conversion.cvcf", std::ios::binary);
    vc::cvcf::reader compact_input(ifs);
    vc::cvcf::marker buff;
    vc::cvcf::reader::input_iterator cur(compact_input, buff);
    vc::cvcf::reader::input_iterator end;

    while (cur != end)
    {
      std::cout << "22" << "\t" << cur->pos() << "\t" << cur->ref() << "\t" << cur->alt() << "\t";

      for (auto gt = cur->begin(); gt != cur->end(); )
      {
        if (*gt == vc::allele_status::has_ref)
          std::cout << "0";
        else if (*gt == vc::allele_status::has_alt)
          std::cout << "1";
        else
          std::cout << ".";

        ++gt;
        if (gt != cur->end())
          std::cout << "\t";
      }

      std::cout << std::endl;

      ++cur;
    }
  }

}

int main(int argc, char** argv)
{
  convert_file_test();
  varint_test();

  return 0;
}


