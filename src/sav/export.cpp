#include <cmath>
#include "sav/export.hpp"
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"
#include "savvy/utility.hpp"

#include <fstream>
#include <ctime>
#include <getopt.h>

class export_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string input_path_;
  std::string output_path_;
  bool help_ = false;
  savvy::fmt format_ = savvy::fmt::allele;
public:
  export_prog_args() :
    long_options_(
      {
        {"format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() { return input_path_; }
  const std::string& output_path() { return output_path_; }
  savvy::fmt format() const { return format_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav export [opts] [in.sav] [out.{vcf,vcf.gz,bcf}]\n";
    os << "\n";
    //os << " -f, --format     : Format field to copy (GT or GP, default: GT)\n";
    os << " -h, --help       : Print usage\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "f:h", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt)
      {
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
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      input_path_ = "/dev/stdin";
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      output_path_ = "/dev/stdout";
    }
    else if (remaining_arg_count == 2)
    {
      input_path_ = argv[optind];
      output_path_ = argv[optind + 1];
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int export_main(int argc, char** argv)
{
  export_prog_args args;
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

  savvy::sav::reader<1> input(args.input_path(), args.format());
  savvy::site_info variant;
  std::vector<float> genotypes;

  std::vector<std::string> sample_ids(input.samples_end() - input.samples_begin());
  std::copy(input.samples_begin(), input.samples_end(), sample_ids.begin());
  auto variant_metadata = input.prop_fields();

  auto headers = input.headers();

  for (auto it = headers.begin(); it != headers.end(); )
  {
    std::string header_id = savvy::parse_header_id(it->second);
    if (it->first == "INFO" && (header_id == "ID" || header_id == "QUAL" || header_id == "FILTER"))
      it = headers.erase(it);
    else
    {
      if (it->first == "FORMAT")
      {
        if (header_id == "GT")
          it->second = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
        else if (header_id == "HDS")
          it->second = "<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage\">";
        else if (header_id == "GP")
          it->second = "<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1\">"; // TODO: Handle other ploidy levels.
      }
      else if (it->first == "fileDate")
      {
        std::time_t t = std::time(nullptr);
        char datestr[11];
        if (std::strftime(datestr, sizeof(datestr), "%Y%m%d", std::localtime(&t)))
        {
          it->second = std::string(datestr);
        }
      }

      ++it;
    }
  }

  savvy::vcf::writer vcf_output(args.output_path(), sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end());

  while (input.read(variant, genotypes))
    vcf_output.write(variant, genotypes);

  return 0;
}