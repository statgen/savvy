/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "sav/export.hpp"
#include "sav/sort.hpp"
#include "sav/utility.hpp"
#include "sav/filter.hpp"

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
  std::vector<std::string> info_fields_;
  filter filter_;
  std::string input_path_;
  std::string output_path_;
  std::string index_path_;
  std::string file_format_;
  std::string headers_path_;
  std::unique_ptr<savvy::s1r::sort_point> sort_type_;
  savvy::fmt format_ = savvy::fmt::gt;
  savvy::bounding_point bounding_point_ = savvy::bounding_point::beg;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  bool help_ = false;
  bool index_ = false;
public:
  export_prog_args() :
    long_options_(
      {
        {"block-size", required_argument, 0, 'b'},
        {"bounding-point", required_argument, 0, 'p'},
        {"data-format", required_argument, 0, 'd'},
        {"file-format", required_argument, 0, 'f'},
        {"filter", required_argument, 0, 'e'},
        {"headers", required_argument, 0, '\x01'},
        {"help", no_argument, 0, 'h'},
        {"index", no_argument, 0, 'x'},
        {"index-file", required_argument, 0, 'X'},
        {"info-fields", required_argument, 0, 'm'},
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
  const std::string& index_path() const { return index_path_; }
  const std::string& headers_path() const { return headers_path_; }
  const filter& filter_functor() const { return filter_; }
  const std::set<std::string>& subset_ids() const { return subset_ids_; }
  const std::vector<savvy::region>& regions() const { return regions_; }
  const std::vector<std::string>& info_fields() const { return info_fields_; }
  const std::unique_ptr<savvy::s1r::sort_point>& sort_type() const { return sort_type_; }
  savvy::fmt format() const { return format_; }
  savvy::bounding_point bounding_point() const { return bounding_point_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  bool index_is_set() const { return index_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav export [opts ...] [in.sav] [out.{vcf,vcf.gz,sav}]\n";
    os << "\n";
    os << " -#                     Number (#) of compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size       Number of markers in SAV compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -d, --data-format      Format field to export (GT, DS, HDS or GP, default: GT)\n";
    os << " -e, --filter           Expression for filtering based on info fields (eg, -e 'AC>=10;AF>0.01') # (IN DEVELOPMENT) More complex expressions in the works\n";
    os << " -f, --file-format      File format (vcf, vcf.gz or sav, default: vcf)\n";
    os << " -h, --help             Print usage\n";
    os << " -i, --sample-ids       Comma separated list of sample IDs to subset\n";
    os << " -I, --sample-ids-file  Path to file containing list of sample IDs to subset\n";
    os << " -m, --info-fields      Comma separated list of INFO (metadata) fields to include with each variant (default: exports all info fields)\n";
    os << " -p, --bounding-point   Determines the inclusion policy of indels during region queries (any, all, beg or end, default: beg)\n";
    os << " -r, --regions          Comma separated list of regions formatted as chr[:start-end]\n";
    os << " -R, --regions-file     Path to file containing list of regions formatted as chr<tab>start<tab>end\n";
    os << " -s, --sort             Enables sorting by first position of allele\n";
    os << " -S, --sort-point       Enables sorting and specifies which allele position to sort by (beg, mid or end)\n";
    os << " -x, --index            Enables indexing (SAV output only)\n";
    os << " -X, --index-file       Enables indexing and specifies index output file (SAV output only)\n";
    os << "\n";
    os << "     --headers          Path to headers file that is either formated as VCF headers or tab-delimited key value pairs\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:d:e:f:hi:I:m:p:r:R:sS:xX:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '\x01':
        {
          if (std::string(long_options_[long_index].name) == "headers")
          {
            headers_path_ = std::string(optarg ? optarg : "");
            break;
          }
          std::cerr << "Invalid long only index (" << long_index << ")\n";
          return false;
        }
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (compression_level_ < 0)
            compression_level_ = 0;
          compression_level_ *= 10;
          compression_level_ += copt - '0';
          break;
        case 'b':
          block_size_ = std::uint16_t(std::atoi(optarg) > 0xFFFF ? 0xFFFF : std::atoi(optarg));
          break;
        case 'd':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          if (str_opt_arg == "HDS")
          {
            format_ = savvy::fmt::hds;
          }
          else if (str_opt_arg == "DS")
          {
            format_ = savvy::fmt::ds;
          }
          else if (str_opt_arg == "GP")
          {
            format_ = savvy::fmt::gp;
          }
          else if (str_opt_arg != "GT")
          {
            std::cerr << "Invalid format field value (" << str_opt_arg << ")\n";
            return false;
          }
          break;
        }
        case 'e':
        {
          std::string str_opt_arg(optarg ? optarg : "");
          filter_ = str_opt_arg;
          if (!filter_)
          {
            std::cerr << "Invalid filter expression (" << str_opt_arg << ")\n";
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
          return true;
        case 'i':
          subset_ids_ = split_string_to_set(optarg ? optarg : "", ',');
          break;
        case 'I':
          subset_ids_ = split_file_to_set(optarg ? optarg : "");
          break;
        case 'm':
          info_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
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
        case 'x':
          index_ = true;
          break;
        case 'X':
          index_ = true;
          index_path_ = optarg;
        default:
          return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 2 && index_ && index_path_.empty())
    {
      std::cerr << "--index-file must be specified when output path is not." << std::endl;
      return false;
    }

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

      if (index_ && index_path_.empty())
        index_path_ = output_path_ + ".s1r";

      if (::savvy::detail::has_extension(output_path_, ".sav"))
        file_format_ = "sav";
      else if (::savvy::detail::has_extension(output_path_, ".vcf"))
        file_format_ = "vcf";
      else if (::savvy::detail::has_extension(output_path_, ".vcf.gz"))
        file_format_ = "vcf.gz";
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (compression_level_ < 0)
      compression_level_ = default_compression_level;
    else if (compression_level_ > 19)
      compression_level_ = 19;

    if (info_fields_.size() && file_format_ == "sav")
    {
      info_fields_.reserve(info_fields_.size() + 3);
      info_fields_.emplace_back("ID");
      info_fields_.emplace_back("QUAL");
      info_fields_.emplace_back("FILTER");
    }

    return true;
  }
};

//template <typename Rdr>
//std::function<bool(const savvy::site_info& var)> gen_filter_predicate(const Rdr& in, const export_prog_args& args)
//{
//  std::size_t hap_cnt = in.samples().size() * in.ploidy();
//  if (args.exclude_monomorphic())
//  {
//    return [hap_cnt](const savvy::site_info& var)
//    {
//      if (var.prop("AC").size())
//      {
//        std::size_t ac = (std::size_t) std::atol(var.prop("AC").c_str());
//        if (!ac || ac == hap_cnt)
//          return false;
//      }
//      return true;
//    };
//  }
//  else
//  {
//    return [](const savvy::site_info& var) { return true; };
//  }
//}

template <typename Vec, typename Writer>
int export_records(savvy::sav::reader& in, const std::vector<savvy::region>& regions, const filter& fn, Writer& out)
{
  savvy::site_info variant;
  Vec genotypes;


  //auto fn = gen_filter_predicate(in, args);

  while (in.read_if(std::ref(fn), variant, genotypes))
    out.write(variant, genotypes);

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Vec, typename Writer>
int export_records(savvy::sav::indexed_reader& in, const std::vector<savvy::region>& regions, const filter& fn, Writer& out)
{
  savvy::site_info variant;
  Vec genotypes;

  //auto fn = gen_filter_predicate(in, args);

  while (in.read_if(std::ref(fn), variant, genotypes))
    out.write(variant, genotypes);

  if (regions.size())
  {
    for (auto it = regions.begin() + 1; it != regions.end(); ++it)
    {
      in.reset_region(*it);
      while (in.read_if(std::ref(fn), variant, genotypes))
        out.write(variant, genotypes);
    }
  }

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Vec, typename Rdr, typename Wrtr>
int prep_writer_for_export(Rdr& input, Wrtr& output, const std::vector<std::string>& sample_ids, const std::vector<std::pair<std::string, std::string>>& headers, const export_prog_args& args)
{
  if (args.sort_type())
  {
    return (sort_and_write_records<std::vector<float>>((*args.sort_type()), input, input.data_format(), args.regions(), output, args.format()) && !input.bad() ? EXIT_SUCCESS : EXIT_FAILURE);
  }
  else
  {
    return export_records<Vec>(input, args.regions(), args.filter_functor(), output);
  }
}

template <typename T>
int prep_reader_for_export(T& input, const export_prog_args& args)
{
  if (!input)
  {
    std::cerr << "Could not open file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }

  std::vector<std::string> sample_ids(input.samples().size());
  std::copy(input.samples().begin(), input.samples().end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = input.subset_samples(args.subset_ids());

  std::vector<std::pair<std::string, std::string>> headers;
  if (args.headers_path().empty())
    headers = input.headers();
  else
  {
    std::ifstream headers_ifs(args.headers_path(), std::ios::binary);
    if (!headers_ifs)
    {
      std::cerr << "Could not open headers file (" << args.headers_path() << ")\n";
      return EXIT_FAILURE;
    }

    std::string line;
    while (std::getline(headers_ifs, line))
    {
      if (line.size())
      {
        line.erase(0, line.find_first_not_of('#'));
        auto delim = line.find_first_of("=\t");
        std::string key = line.substr(0, delim);
        std::string value;
        if (delim != std::string::npos)
          value = line.substr(delim + 1);
        if (std::min(key.size(), value.size()) == 0)
        {
          std::cerr << "Invalid header in " << args.headers_path() << "\n";
          return EXIT_FAILURE;
        }
        headers.emplace_back(std::move(key), std::move(value));
        headers.reserve((std::size_t) headers.size() * 1.5f);
      }
    }
  }


  for (auto it = headers.begin(); it != headers.end(); )
  {
    std::string header_id = savvy::parse_header_sub_field(it->second, "ID");
    if ((it->first == "FORMAT") ||
      (it->first == "INFO"  && args.file_format() != "sav" && (header_id == "ID" || header_id == "QUAL" || header_id == "FILTER")) ||
      (it->first == "INFO" && args.info_fields().size() && std::find(args.info_fields().begin(), args.info_fields().end(), header_id) == args.info_fields().end()))
    {
      it = headers.erase(it);
    }
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
    savvy::sav::writer::options opts;
    opts.compression_level = args.compression_level();
    opts.block_size = args.block_size();
    if (args.index_path().size())
      opts.index_path = args.index_path();

    savvy::sav::writer output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());
    return prep_writer_for_export<savvy::compressed_vector<float>>(input, output, sample_ids, headers, args);
  }
  else
  {
    savvy::vcf::writer<1>::options opts;
    if (args.file_format() == "vcf.gz")
      opts.compression = savvy::vcf::compression_type::bgzip;
    savvy::vcf::writer<1> output(args.output_path(), opts, sample_ids.begin(), sample_ids.end(), headers.begin(), headers.end(), args.format());
    return prep_writer_for_export<std::vector<float>>(input, output, sample_ids, headers, args);
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