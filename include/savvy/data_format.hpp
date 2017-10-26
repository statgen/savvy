
#ifndef LIBSAVVY_DATA_FORMAT_HPP
#define LIBSAVVY_DATA_FORMAT_HPP

#include <cstdint>

namespace savvy
{
  enum class fmt : std::uint8_t
  {
    allele = 1,
    genotype,
    genotype_probability,
    genotype_likelihood,
    phred_scaled_genotype_likelihood,
    dosage,
    haplotype_dosage,
    //phase,
//    gt = genotype,
//    gp = genotype_probability,
//    gl = genotype_likelihood,
//    pl = phred_scaled_genotype_likelihood,
//    ds = dosage
//    ec = dosage
  };
}
#endif //LIBSAVVY_DATA_FORMAT_HPP
