/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

//
//#include "sav_reader.hpp"
//
//#include <fstream>
//#include <sav_reader.hpp>
//#include <allele_vector.hpp>
//
//struct sorting_info
//{
//  std::uint64_t sort_weight;
//  std::size_t original_index;
//
//  sorting_info(std::size_t index) : sort_weight(0), original_index(index) {}
//};
//
int main(int argc, char** argv)
{
//  const std::uint8_t sort_weight_coefficient = 0;
//
//  savvy::sav::reader input(argv[1]);
//
//  const std::size_t sample_count = input.sample_count();
//  const std::size_t ploidy_level = input.ploidy();
//  const std::size_t haplotype_count = sample_count * ploidy_level;
//
//  std::vector<sorting_info> alt_cnts;
//  alt_cnts.reserve(sample_count);
//  for (std::size_t i = 0; i < sample_count; ++i)
//    alt_cnts.emplace_back(i);
//
//  std::for_each(savvy::sav::dense_variant_iterator<float>(input), savvy::sav::dense_variant_iterator<float>(), [&alt_cnts, ploidy_level](const savvy::dense_allele_vector<float>& mkr)
//  {
//
//    const double af = mkr.calculate_allele_frequency();
//    for (auto it = mkr.non_ref_begin(); it != mkr.non_ref_end(); ++it)
//    {
//      alt_cnts[it->offset / ploidy_level].sort_weight += (sort_weight_coefficient * af) + 1;
//    }
//
//  });
//
//  std::sort(alt_cnts.begin(), alt_cnts.end(), [](const sorting_info& a, const sorting_info& b)
//  {
//    return a.sort_weight > b.sort_weight;
//  });
//
//  std::vector<std::size_t> old_order_to_new_order_mapping(sample_count);
//  for (std::size_t i = 0; i < sample_count; ++i)
//    old_order_to_new_order_mapping[alt_cnts[i].original_index] = i;
//
//  std::vector<std::string> sample_ids(sample_count);
//
//  {
//    std::size_t i = 0;
//    for (auto it = input.samples_begin(); it != input.samples_end(); ++it,++i)
//      sample_ids[old_order_to_new_order_mapping[i]] = *it;
//
//  }
//
//  savvy::sav::writer compact_output(argv[2], input.chromosome(), input.ploidy(), sample_ids.begin(), sample_ids.end());
//  savvy::sav::reader input2(argv[1]);
//  std::for_each(savvy::sav::reader::input_iterator(input2, buff), savvy::sav::reader::input_iterator(), [&old_order_to_new_order_mapping, &compact_output, haplotype_count, ploidy_level](const savvy::sav::marker& mrkr)
//  {
//    auto gt_beg = mrkr.non_ref_begin();
//    auto gt_end = mrkr.non_ref_end();
//
//    std::vector<savvy::sav::marker::sparse_vector_allele> sparse_alleles;
//    sparse_alleles.reserve(gt_end - gt_beg);
//
//    for (std::size_t i = 0; gt_beg != gt_end; ++i,++gt_beg)
//    {
//      savvy::sav::marker::sparse_vector_allele tmp(gt_beg->status, (old_order_to_new_order_mapping[gt_beg->offset / ploidy_level] * ploidy_level) + (gt_beg->offset % ploidy_level));
//
//      auto find_res = std::upper_bound( sparse_alleles.begin(), sparse_alleles.end(), tmp, [](const savvy::sav::marker::sparse_vector_allele& a, const savvy::sav::marker::sparse_vector_allele& b)
//      {
//        return (a.offset < b.offset);
//      });
//
//      sparse_alleles.emplace(find_res, std::move(tmp));
//    }
//
//    compact_output << savvy::sav::marker(mrkr.pos(), mrkr.ref(), mrkr.alt(), sparse_alleles.begin(), sparse_alleles.end(), haplotype_count);
//  });
//
  return 0;
}