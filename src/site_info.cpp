#include "savvy/site_info.hpp"

namespace savvy
{
  const std::string site_info::empty_string = {};
  namespace detail
  {
    void print_vcf_site_info(std::ostream& out, const site_info& in, const std::vector<std::string>& info_fields)
    {
      out << in.chromosome() << "\t"
          << in.locus() << "\t"
          << (in.prop("ID").size() ? in.prop("ID") : ".") << "\t"
          << in.ref() << "\t"
          << in.alt() << "\t"
          << (in.prop("QUAL").size() ? in.prop("QUAL") : ".") << "\t"
          << (in.prop("FILTER").size() ? in.prop("FILTER") : ".") << "\t";

      std::size_t counter = 0;
      for (auto it = info_fields.begin(); it != info_fields.end(); ++it)
      {
        std::string val = in.prop(*it);
        if (val.size())
        {
          out << (counter > 0 ? ";" : "") << *it << "=" << val;
          ++counter;
        }
      }

      out << (counter == 0 ? ".\t" : "\t");

      out << "GT";
    }
  }
}