#include "sav/import.hpp"
#include "sav/export.hpp"
#include "savvy/utility.hpp"

#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <iostream>

class prog_args
{
private:
  static const int default_compression_level = 3;
  static const int default_block_size = 2048;

  std::vector<option> long_options_;
  std::string sub_command_;
  bool help_ = false;
  bool version_ = false;
public:
  prog_args() :
    long_options_(
      {
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
      })
  {
  }

  const std::string& sub_command() { return sub_command_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  void print_usage(std::ostream& os)
  {
    os << "----------------------------------------------\n";
    os << "Usage: sav [args] [sub-command]\n";
    os << "\n";
    os << " -h, --help       : Print usage\n";
    os << " -v, --version    : Print version\n";
    os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int& argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, "hv", long_options_.data(), &long_index )) != -1)
    {
      std::string str_opt_arg(optarg ? optarg : "");
      char copt = char(opt & 0xFF);
      switch (copt)
      {
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

    if (remaining_arg_count > 0)
    {
      char* tmp = argv[optind];
      sub_command_ = tmp;
      int i = optind;
      for ( ; (i + 1) < argc; ++i)
      {
        argv[i] = argv[i + 1];
      }
      argv[i] = tmp;
      --argc;
    }
    else if (!help_is_set() && !version_is_set())
    {
      std::cerr << "Missing sub-command\n";
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

  optind = 1;

  if (args.sub_command() == "import")
  {
    return import_main(argc, argv);
  }
  else if (args.sub_command() == "export")
  {
    return export_main(argc, argv);
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "vcf2sav v" << savvy::savvy_version() << std::endl;
    return EXIT_SUCCESS;
  }

  std::cerr << "Invalid sub-command (" << args.sub_command() << ")" << std::endl;
  args.print_usage(std::cerr);

  return EXIT_FAILURE;
}