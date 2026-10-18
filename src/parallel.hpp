#pragma once
#include <algorithm>
#include <atomic>
#include <future>
#include <stdexcept>
#include <vector>

namespace leptonff::detail {
template<class F> void parallel_for(std::size_t size,unsigned threads,const F& f) {
  if(!threads || threads>1024)throw std::invalid_argument("Worker count must be in [1,1024]");
  if(threads==1 || size<2) {for(std::size_t i=0;i<size;++i)f(i);return;}
  const auto count=std::min<std::size_t>(threads,size);
  std::atomic<std::size_t> next{0};std::atomic<bool> stop{false};
  std::vector<std::future<void>> jobs;
  for(std::size_t k=0;k<count;++k)jobs.push_back(std::async(std::launch::async,[&]{
    try {while(!stop) {const auto i=next.fetch_add(1);if(i>=size)break;f(i);}}
    catch(...) {stop=true;throw;}
  }));
  std::exception_ptr error;
  for(auto& job:jobs)try {job.get();}catch(...) {if(!error)error=std::current_exception();}
  if(error)std::rethrow_exception(error);
}
} // namespace leptonff::detail
