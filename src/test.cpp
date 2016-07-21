#include <iostream>
#include <fstream>
#include <algorithm>
#include "cvcf_reader.hpp"
#include "m3vcf_reader.hpp"
#include "vcf_reader.hpp"
#include "test_class.hpp"

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

    std::for_each(it->begin(), it->end(), [](const typename T::input_iterator::value_type::const_iterator::value_type& s)
    {
      vc::allele_status foo = s;
    });

    ++it;
  }

}

int main(int argc, char** argv)
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
      vc::vcf::marker m;
      vc::vcf::reader input(file_path);
      vc::vcf::reader::input_iterator eof;
      vc::vcf::reader::input_iterator cur(input, m);

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

  std::for_each(m.non_ref_begin(), m.non_ref_end(), [&zero_one_two_vec, ploidy_level](const vc::cvcf::marker::sparse_allele& a)
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


