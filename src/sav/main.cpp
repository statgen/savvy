/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "sav/concat.hpp"
#include "sav/export.hpp"
#include "sav/head.hpp"
#include "sav/import.hpp"
#include "sav/index.hpp"
#include "sav/merge.hpp"
#include "sav/rehead.hpp"
#include "sav/sort.hpp"
#include "sav/stat.hpp"
#include "savvy/utility.hpp"

#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>

class prog_args
{
private:
  std::string sub_command_;
  bool help_ = false;
  bool version_ = false;
public:
  prog_args()
  {
  }

  const std::string& sub_command() const { return sub_command_; }
  bool help_is_set() const { return help_; }
  bool version_is_set() const { return version_; }

  void print_usage(std::ostream& os)
  {
    //os << "----------------------------------------------\n";
    os << "Usage: sav <sub-command> [args ...]\n";
    os << "or: sav [opts ...]\n";
    os << "\n";
    os << "Sub-commands:\n";
    os << " export:      Exports SAV to VCF or SAV\n";
    os << " head:        Prints SAV headers or samples IDs\n";
    os << " import:      Imports VCF or BCF into SAV\n";
    os << " index:       Indexes SAV file\n";
    os << " merge:       Merges multiple files into one\n";
    os << " rehead:      Replaces headers without recompressing variant blocks.\n";
    os << " stat-index:  Gathers statistics on s1r index\n";
    os << "\n";
    os << "Options:\n";
    os << " -h, --help     Print usage\n";
    os << " -v, --version  Print version\n";
    //os << "----------------------------------------------\n";
    os << std::flush;
  }

  bool parse(int& argc, char**& argv)
  {
    if (argc > 1)
    {
      std::string str_opt_arg(argv[1]);

      if (str_opt_arg == "-h" || str_opt_arg == "--help")
        help_ = true;
      else if (str_opt_arg == "-v" || str_opt_arg == "--version")
        version_ = true;
      else if (str_opt_arg.size() && str_opt_arg.front() != '-')
      {
        sub_command_ = str_opt_arg;
        --argc;
        ++argv;
      }
    }

    if (!help_is_set() && !version_is_set() && sub_command_.empty())
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

  if (args.sub_command() == "concat")
  {
    return concat_main(argc, argv);
  }
  else if (args.sub_command() == "export")
  {
    return export_main(argc, argv);
  }
  else if (args.sub_command() == "head")
  {
    return head_main(argc, argv);
  }
  else if (args.sub_command() == "import")
  {
    return import_main(argc, argv);
  }
  else if (args.sub_command() == "index")
  {
    return index_main(argc, argv);
  }
  else if (args.sub_command() == "merge")
  {
    return merge_main(argc, argv);
  }
  else if (args.sub_command() == "rehead")
  {
    return rehead_main(argc, argv);
  }
  else if (args.sub_command() == "stat-index")
  {
    return stat_index_main(argc, argv);
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "sav v" << savvy::savvy_version() << std::endl;
    return EXIT_SUCCESS;
  }

  std::cerr << "Invalid sub-command (" << args.sub_command() << ")" << std::endl;
  args.print_usage(std::cerr);

  return EXIT_FAILURE;
}