#include "api_error.hpp"
#include "leptonff/fragmentation_c.h"
#include "leptonff/version.hpp"
namespace leptonff::detail { thread_local std::string last_error; }
namespace leptonff { const char* version() noexcept { return LEPTONFF_VERSION; } }
extern "C" const char* leptonff_last_error(void) {
  return leptonff::detail::last_error.c_str();
}
extern "C" const char* leptonff_version(void) { return leptonff::version(); }
