
#include "vc.hpp"
namespace savvy
{
  namespace detail
  {
    bool has_extension(const std::string& fullString, const std::string& ext)
    {
      if (fullString.length() >= ext.length())
        return (0 == fullString.compare (fullString.length() - ext.length(), ext.length(), ext));
      else
        return false;
    }
  }
}