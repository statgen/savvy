/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "sav/export.hpp"
#include "sav/sort.hpp"
#include "sav/utility.hpp"
#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"

#include <set>
#include <fstream>
#include <ctime>
#include <getopt.h>

class export_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;


  std::vector<option> long_options_;
  std::set<std::string> subset_ids_;
  std::vector<savvy::region> regions_;
  std::string input_path_;
  std::string output_path_;
  std::string file_format_;
  std::unique_ptr<savvy::s1r::sort_point> sort_type_;
  bool help_ = false;
  savvy::fmt format_ = savvy::fmt::allele;
  savvy::bounding_point bounding_point_ = savvy::bounding_point::beg;
public:
  export_prog_args() :
    long_options_(
      {
        {"bounding-point", required_argument, 0, 'p'},
        {"data-format", required_argument, 0, 'd'},
        {"file-format", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"regions", required_argument, 0, 'r'},
        {"regions-file", required_argument, 0, 'R'},
        {"sample-ids", required_argument, 0, 'i'},
        {"sample-ids-file", required_argument, 0, 'I'},
        {"sort", no_argument, 0, 's'},
        {"sort-point", required_argument, 0, 'S'},
        {0, 0, 0, 0}
      }),
    file_format_("vcf")
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& file_format() const { return file_format_; }
  const std::set<std::string>& subset_ids() const { return subset_ids_; }
  const std::vector<savvy::region>& regions() const { return regions_; }
  const std::unique_ptr<savvy::s1r::sort_point>& sort_type() const { return sort_type_; }
  savvy::fmt format() const { return format_; }
  savvy::bounding_point bounding_point() const { return bounding_point_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav export [opts ...] [in.sav] [out.{vcf,sav}]\n";
    os << "\n";
    os << " -d, --data-format     : Format field to export (GT, DS, HDS or GP, default: GT)\n";
    //os << " -e, --filter-expression : File format (vcf or sav, default: vcf)\n";
    os << " -f, --file-format     : File format (vcf, vcf.gz or sav, default: vcf)\n";
    os << " -h, --help            : Print usage\n";
    os << " -i, --sample-ids      : Comma separated list of sample IDs to subset\n";
    os << " -I, --sample-ids-file : Path to file containing list of sample IDs to subset\n";
    os << " -p, --bounding-point  : Determines the inclusion policy of indels during region queries (any, all, beg or end, default: beg)\n";
    os << " -r, --regions         : Comma separated list of regions formatted as chr[:start-end]\n";
    os << " -R, --regions-file    : Path to file containing list of regions formatted as chr<tab>start<tab>end\n";
    os << " -s, --sort            : Enables sorting by first position of allele\n";
    os << " -S, --sort-point      : Enables sorting and specifies which allele position to sort by (beg, mid or end)\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "d:f:hi:I:p:r:R:sS:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case 'd':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::haplotype_dosage;
          }
          else if (str_opt_arg == "DS")
          {
            format_ = savvy::fmt::dosage;
          }
          else if (str_opt_arg == "GP")
          {
            format_ = savvy::fmt::genotype_probability;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 'f':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "sav" || str_opt_arg == "vcf" || str_opt_arg == "vcf.gz")
          {
            file_format_ = str_opt_arg;
          }
          else
          {
            std::cerr << "Invalid file format value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 'h':
          help_ = true;
          break;
        case 'i':
          subset_ids_ = split_string_to_set(optarg, ',');
          break;
        case 'I':
          subset_ids_ = split_file_to_set(optarg);
          break;
        case 'p':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "any")
          {
            bounding_point_ = savvy::bounding_point::any;
          }
          else if (str_opt_arg == "all")
          {
            bounding_point_ = savvy::bounding_point::all;
          }
          else if (str_opt_arg == "beg")
          {
            bounding_point_ = savvy::bounding_point::beg;
          }
          else if (str_opt_arg == "end")
          {
            bounding_point_ = savvy::bounding_point::end;
          }
          else
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 'r':
          for (const auto& r : split_string_to_vector(optarg, ','))
            regions_.emplace_back(string_to_region(r));
          break;
        case 'R':
          for (const auto& r : split_file_to_vector(optarg))
          {
            std::string s = r;
            std::size_t pos = s.find('\t');
            if (pos != std::string::npos)
            {
              s[pos] = ':';
              pos = s.find('\t', pos + 1);
              if (pos != std::string::npos)
                s[pos] = '-';
            }
            regions_.emplace_back(string_to_region(s));
          }
          break;
        case 's':
          sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
          break;
        case 'S':
        {
          std::string sort_str(optarg);
          if (sort_str.size())
          {
            if (sort_str.front()=='b')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::beg);
            }
            else if (sort_str.front()=='e')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::end);
            }
            else if (sort_str.front()=='m')
            {
              sort_type_ = savvy::detail::make_unique<savvy::s1r::sort_point>(savvy::s1r::sort_point::mid);
            }
            else
            {
              std::cerr << "Invalid --sort-point argument (" << sort_str << ")." << std::endl;
              return false;
            }
          }
          break;
        }
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 0)
    {
      if (regions_.size())
      {
        std::cerr << "Input path must be specified when using --regions option." << std::endl;
        return false;
      }

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

template <typename Writer>
int export_records(savvy::sav::reader& in, const std::vector<savvy::region>& regions, Writer& out)
{
  savvy::site_info variant;
  std::vector<float> genotypes;

  while (in.read(variant, genotypes))
    out.write(variant, genotypes);

  return out.good() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Writer>
int export_records(savvy::sav::indexed_reader& in, const std::vector<savvy::region>& regions, Writer& out)
{
  savvy::site_info variant;
  std::vector<float> genotypes;

  while (in.read(variant, genotypes))
    out.write(variant, genotypes);

  if (regions.size())
  {
    for (auto it = regions.begin() + 1; it != regions.end(); ++it)
    {
      in.reset_region(*it);
      while (in.read(variant, genotypes))
        out.write(variant, genotypes);
    }
  }

  return out.good() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Rdr, typename Wrtr>
int prep_writer_for_export(Rdr& input, Wrtr& output, const std::vector<std::string>& sample_ids, const std::vector<std::pair<std::string, std::string>>& headers, const export_prog_args& args)
{
  if (args.sort_type())
  {
    return (sort_and_write_records<std::vector<float>>((*args.sort_type()), input, input.data_format(), args.regions(), output, args.format()) ? EXIT_SUCCESS : EXIT_FAILURE);
  }
  else
  {
    return export_records(input, args.regions(), output);
  }
}

template <typename T>
int prep_reader_for_export(T& input, const export_prog_args& args)
{
  std::vector<std::string> sample_ids(input.samples().size());
  std::copy(input.samples().begin(), input.samples().end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = input.subset_samples(args.subset_ids());

  auto variant_metadata = input.info_fields();

  auto headers = input.headers();

  for (auto it = headers.begin(); it != headers.end(); )
  {
    std::string header_id = savvy::parse_header_id(it->second);
    if ((it->first == "INFO" && (header_id == "ID" || header_id == "QUAL" || header_id == "FILTER")) || it->first == "FORMAT")
      it = headers.erase(it);
    else
    {
      if (it->first == "fileDate")
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

//  if (args.format() == savvy::fmt::allele)
//    headers.emplace_back("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
//  else if (args.format() == savvy::fmt::haplotype_dosage)
//    headers.emplace_back("FORMAT", "<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage\">");
//  else if (args.format() == savvy::fmt::dosage)
//    headers.emplace_back("FORMAT", "<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage\">");
////  else if (args.format() == savvy::fmt::genotype_probability)
////    headers.emplace_back("FORMAT", "<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1\">"); // TODO: Handle other ploidy levels.


  if (args.file_format() == "sav")
  {
    savvy::sav::writer output(args.output_path(), sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());
    return prep_writer_for_export(input, output, sample_ids, headers, args);
  }
  else
  {
    savvy::vcf::writer<1>::options opts;
    if (args.file_format() == "vcf.gz")
      opts.compression = savvy::vcf::compression_type::bgzip;
    savvy::vcf::writer<1> output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());
    return prep_writer_for_export(input, output, sample_ids, headers, args);
  }
}

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

  if (args.regions().size())
  {
    savvy::sav::indexed_reader input(args.input_path(), args.regions().front(), args.bounding_point(), args.format());
    return prep_reader_for_export(input, args);
  }
  else
  {
    savvy::sav::reader input(args.input_path(), args.format());
    return prep_reader_for_export(input, args);
  }
}