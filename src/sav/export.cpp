/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#include "sav/export.hpp"
#include "sav/sort.hpp"
#include "sav/utility.hpp"
#include "sav/filter.hpp"

#include "savvy/vcf_reader.hpp"
#include "savvy/sav_reader.hpp"
#include "savvy/savvy.hpp"
#include "savvy/writer.hpp"
#include "savvy/reader.hpp"

#include <cmath>
#include <set>
#include <fstream>
#include <ctime>
#include <getopt.h>
#include <sys/stat.h>

class export_prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;


  std::vector<option> long_options_;
  std::unordered_set<std::string> subset_ids_;
  std::unordered_set<std::string> fields_to_generate_;
  std::vector<savvy::genomic_region> regions_;
  std::unique_ptr<savvy::slice_bounds> slice_;
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
  int update_info_ = -1;
  int compression_level_ = -1;
  std::uint16_t block_size_ = default_block_size;
  bool sites_only_ = false;
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
        {"generate-info", required_argument, 0, '\x01'},
        {"headers", required_argument, 0, '\x01'},
        {"help", no_argument, 0, 'h'},
        {"index", no_argument, 0, 'x'},
        {"index-file", required_argument, 0, 'X'},
        {"info-fields", required_argument, 0, 'm'},
        {"regions", required_argument, 0, 'r'},
        {"regions-file", required_argument, 0, 'R'},
        {"sample-ids", required_argument, 0, 'i'},
        {"sample-ids-file", required_argument, 0, 'I'},
        {"slice", required_argument, 0, 'c'},
        {"sort", no_argument, 0, 's'},
        {"sort-point", required_argument, 0, 'S'},
        {"sites-only", no_argument, 0, '\x02'},
        {"update-info", required_argument, 0, '\x01'},
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
  const std::unordered_set<std::string>& subset_ids() const { return subset_ids_; }
  const std::unordered_set<std::string>& fields_to_generate() const { return fields_to_generate_; }
  const std::vector<savvy::genomic_region>& regions() const { return regions_; }
  const std::vector<std::string>& info_fields() const { return info_fields_; }
  const std::unique_ptr<savvy::s1r::sort_point>& sort_type() const { return sort_type_; }
  const std::unique_ptr<savvy::slice_bounds>& slice() const { return slice_; }
  savvy::fmt format() const { return format_; }
  savvy::bounding_point bounding_point() const { return bounding_point_; }
  std::uint8_t compression_level() const { return std::uint8_t(compression_level_); }
  std::uint16_t block_size() const { return block_size_; }
  bool update_info() const { return update_info_ != 0; }
  bool index_is_set() const { return index_; }
  bool sites_only_is_set() const { return sites_only_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav export [opts ...] [in.sav] [out.{vcf,vcf.gz,sav}]\n";
    os << "\n";
    os << " -#                     Number (#) of compression level (1-19, default: " << default_compression_level << ")\n";
    os << " -b, --block-size       Number of markers in SAV compression block (0-65535, default: " << default_block_size << ")\n";
    os << " -c, --slice            Range formatted as begin:end (non-inclusive end) that specifies a subset of record offsets within file\n";
    os << " -d, --data-format      Format field to export (GT, DS, HDS or GP, default: GT)\n";
    os << " -e, --filter           Expression for filtering based on info fields (eg, -e 'AC>=10;AF>0.01') # (IN DEVELOPMENT) More complex expressions in the works\n";
    os << " -f, --file-format      File format (vcf, vcf.gz or sav, default: vcf)\n";
    os << " -g, --generate-info    Generate info fields specified as a comma separated list (AC,AN,AF,MAF,SPARSE_OFFSETS_<FMT>,SPARSE_VALUES_<FMT>)\n";
    os << " -h, --help             Print usage\n";
    os << " -i, --sample-ids       Comma separated list of sample IDs to subset\n";
    os << " -I, --sample-ids-file  Path to file containing list of sample IDs to subset\n";
    os << " -m, --info-fields      Comma separated list of INFO (metadata) fields to include with each variant (default: exports all info fields)\n";
    os << " -p, --bounding-point   Determines the inclusion policy of indels during region queries (any, all, beg or end, default: beg)\n";
    os << " -r, --regions          Comma separated list of genomic regions formatted as chr[:start-end]\n";
    os << " -R, --regions-file     Path to file containing list of regions formatted as chr<tab>start<tab>end\n";
    os << " -s, --sort             Enables sorting by first position of allele\n";
    os << " -S, --sort-point       Enables sorting and specifies which allele position to sort by (beg, mid or end)\n";
    os << " -x, --index            Enables indexing (SAV output only)\n";
    os << " -X, --index-file       Enables indexing and specifies index output file (SAV output only)\n";
    os << "\n";
    os << "     --headers          Path to headers file that is either formatted as VCF headers or tab-delimited key value pairs\n";
    os << "     --sites-only       Exclude individual level data.\n";
    //os << "     --update-info      Specifies whether AC, AN, AF and MAF info fields should be updated (always, never or auto, default: auto)\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "0123456789b:c:d:e:f:hi:I:m:p:r:R:sS:xX:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
        case '\x01':
        {
          if (std::string(long_options_[long_index].name) == "generate-info")
          {
            fields_to_generate_ = split_string_to_set(optarg ? optarg : "", ',');
            break;
          }
          else if (std::string(long_options_[long_index].name) == "headers")
          {
            headers_path_ = std::string(optarg ? optarg : "");
            break;
          }
          else if (std::string(long_options_[long_index].name) == "update-info")
          {
            std::string update_info_string(optarg ? optarg : "");
            if (update_info_string == "always")
              update_info_ = 1;
            else if (update_info_string == "never")
              update_info_ = 0;
            else if (update_info_string == "auto")
              update_info_ = -1;
            else
            {
              std::cerr << "Invalid --update-info value (" << update_info_string << ")\n";
              return false;
            }

            break;
          }
          std::cerr << "Invalid long only index (" << long_index << ")\n";
          return false;
        }
        case '\x02':
        {
          if (std::string(long_options_[long_index].name) == "sites-only")
          {
            sites_only_ = true;
          }
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
        case 'c':
        {
          std::string s(optarg ? optarg : "");
          const std::size_t colon_pos = s.find(':');
          if (colon_pos != std::string::npos)
          {
            std::string sstart = s.substr(0, colon_pos);
            std::string send = s.substr(colon_pos + 1, s.size() - colon_pos - 1);
            std::int64_t istart = std::atoll(sstart.c_str());
            std::int64_t iend = std::atoll(send.c_str());
            if (iend > istart)
            {
              slice_ = savvy::detail::make_unique<savvy::slice_bounds>(istart, iend);
              break;
            }
          }
          std::cerr << "Invalid --slice range.\n";
          return false;
        }
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
          struct stat buf;
          if (stat(optarg ? optarg : "", &buf) != 0)
          {
            std::cerr << "Cannot open --sample-ids-file (" << (optarg ? optarg : "") << ")\n";
            return false;
          }
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

    if (regions_.size() && slice_)
    {
      std::cerr << "--regions cannot be combined with --slice\n";
      return false;
    }

    if (update_info_ < 0)
    {
      update_info_ = subset_ids_.size() ? 1 : 0; // Automatically update info fields if samples are subset.
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

    if (sites_only_ && file_format_ != "vcf" && file_format_ != "vcf.gz")
    {
      std::cerr << "--sites-only is only supported for VCF file format\n";
      return false;
    }

    {
      std::set<std::string> allowed = {"AC", "AN", "AF", "MAF"};
      for (auto it = fields_to_generate_.begin(); it != fields_to_generate_.end(); ++it)
      {
        if (allowed.find(*it) == allowed.end())
        {
          if (*it != "SPARSE_OFFSETS_" + savvy::fmt_to_string(format_) && *it != "SPARSE_VALUES_" + savvy::fmt_to_string(format_))
          {
            std::cerr << "Invalid --generate-info value (" << (*it) << ")\n";
            return false;
          }

          if (format_ != savvy::fmt::gt && format_ != savvy::fmt::hds)
          {
            std::cerr << "Only GT and HDS are supported with --generate-info (" + *it + ")\n";
            return false;
          }
        }
      }
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

template <typename T>
std::string create_sparse_offsets_string(const std::vector<T>& genotypes)
{
  std::stringstream ret;
  //ret.reserve(genotypes.size());

  bool first_entry = true;
  const int max_size = std::numeric_limits<std::size_t>::digits10 + 1;
  char buffer[max_size] = {0};
  const T empty_value = T();
  for (std::size_t i = 0; i < genotypes.size(); ++i)
  {
    if (genotypes[i] != empty_value)
    {
      if (!first_entry)
        ret.put(',');
      //std::snprintf(buffer, max_size, "%lu", i);
      ret << i; //buffer;
      first_entry = false;
    }
  }

  return ret.str();
}

template <typename T>
std::string create_sparse_values_string(const std::vector<T>& genotypes)
{
  std::stringstream ret;

  bool first_entry = true;
  for (std::size_t i = 0; i < genotypes.size(); ++i)
  {
    if (genotypes[i] != T())
    {
      if (!first_entry)
        ret.put(',');
      ret << genotypes[i];
      first_entry = false;
    }
  }

  return ret.str();
}

template <typename T>
std::string create_sparse_offsets_string(const savvy::compressed_vector<T>& genotypes)
{
  std::stringstream ret;

  bool first_entry = true;
  for (auto it = genotypes.begin(); it != genotypes.end(); ++it)
  {
    if (!first_entry)
      ret.put(',');
    ret << it.offset();
    first_entry = false;
  }

  return ret.str();
}

template <typename T>
std::string create_sparse_values_string(const savvy::compressed_vector<T>& genotypes)
{
  std::stringstream ret;

  bool first_entry = true;
  for (auto it = genotypes.begin(); it != genotypes.end(); ++it)
  {
    if (!first_entry)
      ret.put(',');
    ret << (*it);
    first_entry = false;
  }

  return ret.str();
}


template <typename Vec>
void set_info(savvy::site_info& variant, const Vec& genotypes, const std::set<std::string>& fields_to_generate, savvy::fmt format)
{
  std::size_t ac = 0, an = 0;
  float af = -1.f, maf = 0.f;
  for (auto it = fields_to_generate.begin(); it != fields_to_generate.end(); ++it)
  {
    if (*it == "AC")
    {
      if (af < 0.f)
        std::tie(ac, an, af, maf) = savvy::generate_standard_info_fields(genotypes);
      variant.prop("AC", std::to_string(ac));
    }
    else if (*it == "AN")
    {
      if (af < 0.f)
        std::tie(ac, an, af, maf) = savvy::generate_standard_info_fields(genotypes);
      variant.prop("AN", std::to_string(an));
    }
    else if (*it == "AF")
    {
      if (af < 0.f)
        std::tie(ac, an, af, maf) = savvy::generate_standard_info_fields(genotypes);
      variant.prop("AF", std::to_string(af));
    }
    else if (*it == "MAF")
    {
      if (af < 0.f)
        std::tie(ac, an, af, maf) = savvy::generate_standard_info_fields(genotypes);
      variant.prop("MAF", std::to_string(maf));
    }
    else if (*it == "SPARSE_OFFSETS_" + savvy::fmt_to_string(format))
    {
      variant.prop("SPARSE_OFFSETS_" + savvy::fmt_to_string(format), create_sparse_offsets_string(genotypes));
    }
    else if (*it == "SPARSE_VALUES_" + savvy::fmt_to_string(format))
    {
      variant.prop("SPARSE_VALUES_" + savvy::fmt_to_string(format), create_sparse_values_string(genotypes));
    }
  }
}

template <typename Vec, typename Writer>
int export_records(savvy::sav::reader& in, const export_prog_args& args, Writer& out)
{
  savvy::site_info variant;
  Vec genotypes;


  //auto fn = gen_filter_predicate(in, args);

  while (in.read_if(std::ref(args.filter_functor()), variant, genotypes))
  {
    if (args.update_info())
      savvy::update_info_fields(variant, genotypes, args.format());

    set_info(variant, genotypes, {args.fields_to_generate().begin(), args.fields_to_generate().end()}, args.format());

    out.write(variant, args.sites_only_is_set() ? Vec() : genotypes);
  }

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Vec, typename Writer>
int export_records(savvy::sav::indexed_reader& in, const export_prog_args& args, Writer& out)
{
  savvy::site_info variant;
  Vec genotypes;

  //auto fn = gen_filter_predicate(in, args);

  while (in.read_if(std::ref(args.filter_functor()), variant, genotypes))
  {
    if (args.update_info())
      savvy::update_info_fields(variant, genotypes, args.format());

    set_info(variant, genotypes, {args.fields_to_generate().begin(), args.fields_to_generate().end()}, args.format());

    out.write(variant, genotypes);
  }

  if (args.regions().size())
  {
    for (auto it = args.regions().begin() + 1; it != args.regions().end(); ++it)
    {
      in.reset_bounds(*it);
      while (in.read_if(std::ref(args.filter_functor()), variant, genotypes))
      {
        if (args.update_info())
          savvy::update_info_fields(variant, genotypes, args.format());

        set_info(variant, genotypes, {args.fields_to_generate().begin(), args.fields_to_generate().end()}, args.format());

        out.write(variant, genotypes);
      }
    }
  }

  return out.good() && !in.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}

template <typename Vec, typename Rdr, typename Wrtr>
int prep_writer_for_export(Rdr& input, Wrtr& output, const std::vector<std::string>& sample_ids, const std::vector<std::pair<std::string, std::string>>& headers, const export_prog_args& args)
{
  if (args.sort_type())
  {
    return (sort_and_write_records<std::vector<float>, less_than_comparator>((*args.sort_type()), input, input.data_format(), args.regions(), output, args.format(), args.update_info()) && !input.bad() ? EXIT_SUCCESS : EXIT_FAILURE);
  }
  else
  {
    return export_records<Vec>(input, args, output);
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
    sample_ids = input.subset_samples({args.subset_ids().begin(), args.subset_ids().end()});

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
        if (headers.capacity() == headers.size())
          headers.reserve((std::size_t) headers.size() * 1.5f);
        headers.emplace_back(std::move(key), std::move(value));
      }
    }
  }

  std::set<std::string> info_fields_already_included;
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
      else if (it->first == "INFO")
      {
        info_fields_already_included.insert(savvy::parse_header_sub_field(it->second, "ID"));
      }

      ++it;
    }
  }

  for (auto it = args.fields_to_generate().begin(); it != args.fields_to_generate().end(); ++it)
  {
    if (info_fields_already_included.find(*it) == info_fields_already_included.end())
      headers.emplace_back("INFO", "<ID=" + (*it) + ">");
  }

  if (args.sites_only_is_set())
  {
    sample_ids.clear();
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

int export_main_old(int argc, char** argv)
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
  else if (args.slice())
  {
    savvy::sav::indexed_reader input(args.input_path(), {""}, args.bounding_point(), args.format());
    input.reset_bounds(*args.slice());
    return prep_reader_for_export(input, args);
  }
  else
  {
    savvy::sav::reader input(args.input_path(), args.format());
    return prep_reader_for_export(input, args);
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

  savvy::v2::reader rdr(args.input_path());
  if (!rdr)
  {
    std::cerr << "Error: failed to open input file" << std::endl;
    return EXIT_FAILURE;
  }

  if (args.regions().size())
  {
    if (rdr.reset_bounds(args.regions().front()).bad())
    {
      std::cerr << "Error: failed to load index for genomic region query" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (args.slice())
  {
    if (rdr.reset_bounds(*args.slice()).bad())
    {
      std::cerr << "Error: failed to load index for slice query" << std::endl;
      return EXIT_FAILURE;
    }
  }

  auto fmt = savvy::file::format::vcf;
  if (args.file_format() == "sav")
    fmt = savvy::file::format::sav2;
  else if (args.file_format() == "bcf")
    fmt = savvy::file::format::bcf;

  std::vector<std::string> sample_ids(rdr.samples().size());
  std::copy(rdr.samples().begin(), rdr.samples().end(), sample_ids.begin());
  if (args.subset_ids().size())
    sample_ids = rdr.subset_samples({args.subset_ids().begin(), args.subset_ids().end()});

  savvy::v2::writer wrt(args.output_path(), fmt, rdr.headers(), sample_ids);

  savvy::v2::variant r;
  while (rdr.read(r))
  {
    if (true) //args.filter_functor()(r))
    {
      wrt.write(r);
    }
  }

  return wrt.good() && !rdr.bad() ? EXIT_SUCCESS : EXIT_FAILURE;
}