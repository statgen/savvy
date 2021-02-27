/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/stat.hpp"
#include "sav/utility.hpp"
#include "savvy/s1r.hpp"
#include "savvy/savvy.hpp"
#include "savvy/reader.hpp"
#include "savvy/writer.hpp"
#include "sav/filter.hpp"

#include <functional>
#include <getopt.h>
#include <memory>

class stat_prog_args
{
private:
  filter filter_;
  std::vector<option> long_options_;
  std::string input_path_;
  std::string summary_path_ = "/dev/stdout";
  std::string per_ac_path_;
  std::string per_sample_path_;
  std::unique_ptr<savvy::genomic_region> reg_;
  bool help_ = false;
public:
  stat_prog_args() :
    long_options_(
      {
        {"filter", required_argument, 0, 'f'},
        {"help", no_argument, 0, 'h'},
        {"per-ac-out", required_argument, 0, '\x01'},
        {"per-sample-out", required_argument, 0, '\x01'},
        {"region", required_argument, 0, 'r'},
        {"summary-out", required_argument, 0, '\x01'},
        {0, 0, 0, 0}
      })
  {
  }

  const filter& filter_functor() const { return filter_; }
  const std::string& input_path() const { return input_path_; }
  const std::string& summary_path() const { return summary_path_; }
  const std::string& per_ac_path() const { return per_ac_path_; }
  const std::string& per_sample_path() const { return per_sample_path_; }
  const std::unique_ptr<savvy::genomic_region>& reg() const { return reg_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav stat [opts ...] <in.sav> \n";
    os << "\n";
    os << " -h, --help  Print usage\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "\x01:hr:", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case '\x01':
      {
        std::string long_opt_name = long_options_[long_index].name;
        if (long_opt_name == "per-ac-out")
        {
          per_ac_path_ = optarg ? optarg : "";
          break;
        }
        else if (long_opt_name == "per-sample-out")
        {
          per_sample_path_ = optarg ? optarg : "";
          break;
        }

        std::cerr << "Invalid long only index (" << long_index << ")\n";
        return false;
      }
      case 'f':
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
      case 'h':
        help_ = true;
        return true;
      case 'r':
        reg_ = savvy::detail::make_unique<savvy::genomic_region>(string_to_region(optarg ? optarg : ""));
        break;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

struct per_ac_t
{
  std::size_t n_snp = 0;
  std::size_t n_indel = 0;
  std::size_t n_syn = 0;
  std::size_t n_nonsyn = 0;

  static void print_header(std::ostream& os)
  {
    os << "#bin_id\tn_snp\tn_indel\tn_syn\tn_nonsyn\n";
  }

  void print(std::ostream& os, std::size_t bin) const
  {
    os << bin << "\t"
       << n_snp  << "\t"
       << n_indel << "\t"
       << n_syn << "\t"
       << n_nonsyn << "\n";
  }
};

struct per_sample_t
{
  std::string sample_id;
  std::size_t n_het = 0;
  std::size_t n_hom = 0;
  std::size_t n_snp = 0;
  std::size_t n_indel = 0;
  std::size_t n_syn = 0;
  std::size_t n_nonsyn = 0;

  per_sample_t(std::string sid) :
    sample_id(std::move(sid))
  {
  }

  static void print_header(std::ostream& os)
  {
    os << "#sample_id\tn_het\tn_hom\tn_snp\tn_indel\tn_syn\tn_nonsyn\n";
  }

  void print(std::ostream& os) const
  {
    os << sample_id << "\t"
       << n_het << "\t"
       << n_hom << "\t"
       << n_snp  << "\t"
       << n_indel << "\t"
       << n_syn << "\t"
       << n_nonsyn << "\n";
  }
};

int stat_main(int argc, char** argv)
{
  stat_prog_args args;
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

  std::size_t multi_allelic{}, record_cnt{}, variant_cnt{};


  savvy::reader input_file(args.input_path());
  if (!input_file)
  {
    std::cerr << "Error: could not open " << args.input_path() << std::endl;
    return EXIT_FAILURE;
  }

  if (args.reg())
  {
    input_file.reset_bounds(*args.reg());
    if (!input_file)
    {
      std::cerr << "Error: could not load region " << args.reg()->chromosome() << ":" << args.reg()->from() << "-" << args.reg()->to() << std::endl;
      return EXIT_FAILURE;
    }
  }

  savvy::variant rec;
  std::vector<std::int8_t> geno;

  std::vector<per_sample_t> per_sample_stats;
  if (args.per_sample_path().size())
    per_sample_stats.assign(input_file.samples().begin(), input_file.samples().end());

  std::size_t bin_width = 1;
  std::vector<per_ac_t> per_ac_stats;

  std::unordered_set<std::string> synonymous_labels = {
    "start_retained",
    "stop_retained",
    "synonymous"};

  std::unordered_set<std::string> nonsynonymous_labels = {
    "stop_gained",
    "frameshift",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense"};

  while (input_file.read(rec))
  {
    if (!args.filter_functor()(rec)) continue;

    if (rec.alts().size() > 1)
      ++multi_allelic;
    variant_cnt += std::max<std::size_t>(1, rec.alts().size());
    ++record_cnt;

    bool is_snp = rec.ref().size() == 1 && rec.alts()[0].size() == 1;
    bool is_syn = false;
    bool is_nonsyn = false;
    std::string ann;
    if (rec.get_info("ANN", ann))
    {
      std::size_t scnt = 0, nonscnt = 0;
      std::vector<std::string> transcripts = split_string_to_vector(ann.c_str(), ',');
      for (auto it = transcripts.begin(); it != transcripts.end(); ++it)
      {
        std::vector<std::string> fields = split_string_to_vector(it->c_str(), '|');
        if (fields.size() >= 2)
        {
          std::vector<std::string> effects = split_string_to_vector(fields[1].c_str(), '&');
          for (auto jt = effects.begin(); jt != effects.end(); ++jt)
          {
            if (synonymous_labels.find(*jt) != synonymous_labels.end())
              ++scnt;
            else if (nonsynonymous_labels.find(*jt) != nonsynonymous_labels.end())
              ++nonscnt;
          }
        }
      }

      //is_syn = scnt && !nonscnt;
      //is_nonsyn = !scnt && nonscnt;
      is_nonsyn = nonscnt > 0;
      is_syn = scnt && !nonscnt;
    }

    if (args.per_ac_path().size())
    {
      std::int64_t ac,an;
      if (!rec.get_info("AC", ac) || !rec.get_info("AN", an))
      {
        std::cerr << "Error: AC and AN INFO fields are required" << std::endl;
        return EXIT_FAILURE;
      }

      if (ac > an || ac < 0)
      {
        std::cerr << "Error: AC INFO field must be in range of [0, AN]" << std::endl;
        return EXIT_FAILURE;
      }

      if (an / bin_width + 1 > per_ac_stats.size())
        per_ac_stats.resize(an / bin_width + 1);

      auto& s = per_ac_stats[ac / bin_width];

      if (is_snp)
        s.n_snp += 1;
      else
        s.n_indel += 1;

      if (is_syn)
        s.n_syn += 1;
      if (is_nonsyn)
        s.n_nonsyn += 1;
    }

    if (per_sample_stats.size())
    {
      if (!rec.get_format("GT", geno)) continue;

      savvy::stride_reduce(geno, geno.size() / per_sample_stats.size());

      for (std::size_t i = 0; i < geno.size(); ++i)
      {
        int8_t g = geno[i];
        if (g < 0) continue;
        if (g)
        {
          if (is_snp)
            per_sample_stats[i].n_snp += g;
          else
            per_sample_stats[i].n_indel += g;

          if (g == 1)
            ++per_sample_stats[i].n_het;
          else // assuming  g == 2
            ++per_sample_stats[i].n_hom;

          if (is_syn)
            per_sample_stats[i].n_syn += g;
          if (is_nonsyn)
            per_sample_stats[i].n_nonsyn += g;
        }
      }
    }
  }


  std::ofstream summary_out(args.summary_path(), std::ios::binary);
  std::cout << record_cnt << "\t" << variant_cnt << "\t" << multi_allelic << "\n";

  if (args.per_sample_path().size())
  {
    std::ofstream per_sample_out(args.per_sample_path(), std::ios::binary);
    per_sample_t::print_header(per_sample_out);
    for (const per_sample_t& s : per_sample_stats)
    {
      s.print(per_sample_out);
    }
  }

  if (args.per_ac_path().size())
  {
    std::ofstream per_ac_out(args.per_ac_path(), std::ios::binary);
    per_ac_t::print_header(per_ac_out);
    for (std::size_t i = 0; i < per_ac_stats.size(); ++i)
    {
      per_ac_stats[i].print(per_ac_out, i / bin_width);
    }
  }

  return EXIT_SUCCESS;
}

class stat_index_prog_args
{
private:
  std::vector<option> long_options_;
  std::string input_path_;
  bool help_ = false;
public:
  stat_index_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav stat-index [opts ...] <in.sav> \n";
    os << "\n";
    os << " -h, --help  Print usage\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "h", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
      if (savvy::detail::file_exists(input_path_ + ".s1r"))
        input_path_ += ".s1r";
    }
    else if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    return true;
  }
};

int stat_index_main(int argc, char** argv)
{
  stat_index_prog_args args;
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

  std::vector<savvy::s1r::index_statistics> stats = savvy::s1r::stat_index(args.input_path());

  if (stats.empty())
  {
    std::cerr << "Could not open index file (" << args.input_path() << ")\n";
    return EXIT_FAILURE;
  }


  std::cout << "chromosome";
  for (auto it = stats.begin(); it != stats.end(); ++it)
    std::cout << "\t" << it->contig;
  std::cout << std::endl;

  std::cout << "tree height";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->tree_height;
  }
  std::cout << std::endl;

  std::cout << "block count";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->block_count;
  }
  std::cout << std::endl;

  std::cout << "record count";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->record_count;
  }
  std::cout << std::endl;

  std::cout << "min position";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->min_position;
  }
  std::cout << std::endl;

  std::cout << "max position";
  for (auto it = stats.begin(); it != stats.end(); ++it)
  {
    std::cout << "\t" << it->max_position;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

// sav stat-index ./test_file_hard.sav.s1r | grep "^marker count" | cut -f 2- | xargs echo | awk 't=0; {for(i=1;i<=NF;i++) t+=$i; print t}'

class stat_merge_prog_args
{
private:
  std::vector<option> long_options_;
  std::vector<std::string> input_paths_;
  bool help_ = false;
public:
  stat_merge_prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::vector<std::string>& input_paths() const { return input_paths_; }
  bool help_is_set() const { return help_; }

  void print_usage(std::ostream& os)
  {
    os << "Usage: sav stat-merge [opts ...] <in1.tsv> <in2.tsv> [other.tsv ...] \n";
    os << "\n";
    os << " -h, --help  Print usage\n";
    os << std::flush;
  }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "h", long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count < 1)
    {
      std::cerr << "Too few arguments\n";
      return false;
    }
    else
    {
      input_paths_.assign(argv + optind, argv + argc);
    }

    return true;
  }
};

std::int64_t str_to_int(const std::string& str)
{
  return std::atoll(str.c_str());
}

int stat_merge_main(int argc, char** argv)
{
  stat_merge_prog_args args;
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

  std::list<std::ifstream> files;
  for (const std::string& p : args.input_paths())
    files.emplace_back(p);


  std::string line;
  std::size_t cnt = 0;
  for (std::ifstream& f : files) // Skip header line
  {
    if (!f)
    {
      std::cerr << "Error: Failed to open " << args.input_paths()[cnt] << std::endl;
      return EXIT_FAILURE;
    }
    std::getline(f, line);
    ++cnt;
  }
  std::cout << line << std::endl;

  std::vector<std::int64_t> agg_data_vec;
  std::vector<std::int64_t> data_vec;
  std::string line_id;

  while (files.size())
  {
    std::fill(agg_data_vec.begin(), agg_data_vec.end(), 0ll);
    for (auto ft = files.begin(); ft != files.end(); ++ft)
    {
      if (!std::getline(*ft, line))
      {
        if (ft != files.begin())
        {
          std::cerr << "Error: Number of lines does not match" << std::endl;
          return EXIT_FAILURE;
        }
        files.clear();
        break;
      }

      auto fields = split_string_to_vector(line.c_str(), '\t');
      if (fields.empty())
      {
        std::cerr << "Error: Empty line" << std::endl;
        return EXIT_FAILURE;
      }

      if (ft == files.begin())
        line_id = fields[0];

      if (fields[0] != line_id)
      {
        std::cerr << "Error: ID field does not match" << std::endl;
        return EXIT_FAILURE;
      }

      data_vec.resize(fields.size() - 1);
      std::transform(fields.begin() + 1, fields.end(), data_vec.begin(), str_to_int);

      if (agg_data_vec.empty())
        agg_data_vec.resize(data_vec.size(), 0ll);

      if (agg_data_vec.size() != data_vec.size())
      {
        std::cerr << "Error: Mismatch in number of columns" << std::endl;
        return EXIT_FAILURE;
      }

      std::transform(agg_data_vec.begin(), agg_data_vec.end(), data_vec.begin(), agg_data_vec.begin(), std::plus<std::size_t>());
    }

    std::cout << line_id;
    for (std::int64_t d : agg_data_vec)
      std::cout << "\t" << d;
    std::cout << "\n";
  }

  return EXIT_SUCCESS;
}