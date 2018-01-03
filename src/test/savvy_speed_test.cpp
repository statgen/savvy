/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "savvy/reader.hpp"
#include "savvy/savvy.hpp"

#include <stdlib.h>
#include <getopt.h>

#include <fstream>
#include <vector>
#include <chrono>

class prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  //std::string output_path_;
  bool help_ = false;
  bool version_ = false;
  savvy::fmt format_ = savvy::fmt::allele;
public:
  prog_args() :
    long_options_(
      {
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() { return input_path_; }
  //const std::string& output_path() { return output_path_; }
  savvy::fmt format() const { return format_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: savvy-speed-test [args] [in.{vcf,vcf.gz,bcf,sav}]\n";
    os << "\n";
    os << " -f, --format     : Format field to copy (GT or GP, default: GT)\n";
    os << " -h, --help       : Print usage\n";
    os << " -v, --version    : Print version\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "f:hv", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt) {
        case 'f':
          if (str_opt_arg == "GP")
          {
            format_ = savvy::fmt::genotype_probability;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        case 'h':
          help_ = true;
          break;
        case 'v':
          version_ = true;
          break;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_path_ = "/dev/stdin";
      //output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      //output_path_ = "/dev/stdout";
    }
//    else if (remaining_arg_count == 2)
//    {
//      input_path_ = argv[optind];
//      output_path_ = argv[optind + 1];
//    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "Savvy Speed Test v" << savvy::savvy_version() << std::endl;
    return EXIT_SUCCESS;
  }

//  auto* hts_file = bcf_open(args.input_path().c_str(), "r");
//  auto* hts_hdr = bcf_hdr_read(hts_file);
//  bcf1_t* hts_rec = bcf_init1();
//
//  std::vector<float> destination;
//  int * gt = nullptr;
//  int gt_sz = 0;
//
//  auto start = std::chrono::high_resolution_clock::now();
//  std::size_t cnt = 0;
//  while (bcf_read(hts_file, hts_hdr, hts_rec) >= 0)
//  {
//    bcf_unpack(hts_rec, BCF_UN_ALL);
//    bcf_get_genotypes(hts_hdr, hts_rec, &(gt), &(gt_sz));
//    destination.resize(0);
//    destination.resize(gt_sz);
//
//    for (std::size_t j = 1; j < hts_rec->n_allele; ++j)
//    {
//      for (std::size_t i = 0; i < gt_sz; ++i)
//      {
//        if (gt[i] == bcf_gt_missing)
//          destination[i] = std::numeric_limits<float>::quiet_NaN();
//        else
//          destination[i] = (bcf_gt_allele(gt[i]) == j ? 1.0f : 0.0f);
//      }
//      ++cnt;
//    }
//  }
//  long elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count();
//  std::cout << "markers: " << cnt << std::endl;
//  std::cout << "elapsed: " << elapsed << " seconds" << std::endl;
//
//  return 0;

  savvy::reader<1> input(args.input_path(), args.format());

  if (input.good())
  {
    savvy::site_info variant;
    std::vector<float> genotypes;

    auto start = std::chrono::high_resolution_clock::now();
    std::size_t cnt = 0;
    while (input.read(variant, genotypes))
    {
      ++cnt;
    }
    long elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "markers: " << cnt << std::endl;
    std::cout << "elapsed: " << elapsed << "s" << std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}