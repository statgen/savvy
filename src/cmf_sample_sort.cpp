
#include "cmf_reader.hpp"

#include <fstream>

struct sorting_info
{
  std::uint64_t sort_weight;
  std::size_t original_index;

  sorting_info(std::size_t index) : sort_weight(0), original_index(index) {}
};

int main(int argc, char** argv)
{
  const std::uint8_t sort_weight_coefficient = 0;

  std::ofstream ofs(argv[2], std::ios::binary);

  vc::cmf::marker buff;
  std::ifstream ifs(argv[1], std::ios::binary);
  vc::cmf::reader input(ifs);

  const std::size_t sample_count = input.sample_count();
  const std::size_t ploidy_level = input.ploidy();
  const std::size_t haplotype_count = sample_count * ploidy_level;

  std::vector<sorting_info> alt_cnts;
  alt_cnts.reserve(sample_count);
  for (std::size_t i = 0; i < sample_count; ++i)
    alt_cnts.emplace_back(i);

  std::for_each(vc::cmf::reader::input_iterator(input, buff), vc::cmf::reader::input_iterator(), [&alt_cnts, ploidy_level](const vc::cmf::marker& mkr)
  {

    const double af = mkr.calculate_allele_frequency();
    for (auto it = mkr.non_ref_begin(); it != mkr.non_ref_end(); ++it)
    {
      alt_cnts[it->offset / ploidy_level].sort_weight += (sort_weight_coefficient * af) + 1;
    }

  });

  std::sort(alt_cnts.begin(), alt_cnts.end(), [](const sorting_info& a, const sorting_info& b)
  {
    return a.sort_weight > b.sort_weight;
  });

  std::vector<std::size_t> old_order_to_new_order_mapping(sample_count);
  for (std::size_t i = 0; i < sample_count; ++i)
    old_order_to_new_order_mapping[alt_cnts[i].original_index] = i;

  std::vector<std::string> sample_ids(sample_count);

  {
    std::size_t i = 0;
    for (auto it = input.samples_begin(); it != input.samples_end(); ++it,++i)
      sample_ids[old_order_to_new_order_mapping[i]] = *it;

  }

  vc::cmf::writer compact_output(ofs, input.chromosome(), input.ploidy(), sample_ids.begin(), sample_ids.end());
  std::ifstream ifs2(argv[1], std::ios::binary);
  vc::cmf::reader input2(ifs2);
  std::for_each(vc::cmf::reader::input_iterator(input2, buff), vc::cmf::reader::input_iterator(), [&old_order_to_new_order_mapping, &compact_output, haplotype_count, ploidy_level](const vc::cmf::marker& mrkr)
  {
    auto gt_beg = mrkr.non_ref_begin();
    auto gt_end = mrkr.non_ref_end();

    std::vector<vc::cmf::marker::sparse_vector_allele> sparse_alleles;
    sparse_alleles.reserve(gt_end - gt_beg);

    for (std::size_t i = 0; gt_beg != gt_end; ++i,++gt_beg)
    {
      vc::cmf::marker::sparse_vector_allele tmp(gt_beg->status, (old_order_to_new_order_mapping[gt_beg->offset / ploidy_level] * ploidy_level) + (gt_beg->offset % ploidy_level));

      auto find_res = std::upper_bound( sparse_alleles.begin(), sparse_alleles.end(), tmp, [](const vc::cmf::marker::sparse_vector_allele& a, const vc::cmf::marker::sparse_vector_allele& b)
      {
        return (a.offset < b.offset);
      });

      sparse_alleles.emplace(find_res, std::move(tmp));
    }

    compact_output << vc::cmf::marker(mrkr.pos(), mrkr.ref(), mrkr.alt(), sparse_alleles.begin(), sparse_alleles.end(), haplotype_count);
  });

  return 0;
}