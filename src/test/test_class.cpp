/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

//#include "allele_status.hpp"
//#include "cvcf_reader.hpp"
//#include "m3vcf_reader.hpp"
//#include "test_class.hpp"

//template <typename Reader>
//some_analysis<Reader>::some_analysis(Reader& file_reader) :
//  file_reader_(file_reader)
//{
//
//}
//
//template <typename Reader>
//void some_analysis<Reader>::run()
//{
//  typename Reader::input_iterator::buffer buff;
//  typename Reader::input_iterator end;
//  typename Reader::input_iterator it(file_reader_, buff);
//  while (it != end)
//  {
//    double af = it->calculate_allele_frequency();
//    handle_marker(*it);
//    ++it;
//  }
//}
//
//template <typename Reader>
//void some_analysis<Reader>::handle_marker(const typename Reader::input_iterator::value_type& m)
//{
//  for (auto jt = m.begin(); jt != m.end(); ++jt)
//  {
//    savvy::allele_status foo = *jt;
//  }
//}
//
//template <>
//void some_analysis<savvy::sav::reader>::handle_marker(const savvy::sav::marker& m)
//{
//  for (auto jt = m.non_ref_begin(); jt != m.non_ref_end(); ++jt)
//  {
//    std::uint64_t foo = jt->offset;
//    savvy::allele_status bar = jt->status;
//  }
//}
//
//template <typename Reader>
//some_analysis<Reader> make_analysis(Reader& file_reader)
//{
//  some_analysis<Reader> ret(file_reader);
//  return ret;
//}

//template class some_analysis<savvy::sav::reader>;
//template class some_analysis<savvy::m3vcf::reader>;
//template some_analysis<savvy::sav::reader> make_analysis<savvy::sav::reader>(savvy::sav::reader&);
//template some_analysis<savvy::m3vcf::reader> make_analysis<savvy::m3vcf::reader>(savvy::m3vcf::reader&);


