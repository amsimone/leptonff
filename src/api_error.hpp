#pragma once
#include <exception>
#include <string>
namespace leptonff::detail {
extern thread_local std::string last_error;
template<class F> int guarded(F&& action) noexcept {
  try { last_error.clear(); action(); return 0; }
  catch (const std::exception& e) { last_error=e.what(); return 1; }
  catch (...) { last_error="Unknown leptonff failure"; return 1; }
}
} // namespace leptonff::detail
