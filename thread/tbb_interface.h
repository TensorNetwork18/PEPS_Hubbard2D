#ifndef _TBB_INTERFACE_H_
#define _TBB_INTERFACE_H_

#include <memory>
#include <functional>
#include <tbb/tbb.h>


template<typename T = std::function<void()>>
struct compute : public tbb::task {
  private:
    T call;
    task* execute() {
      call();
      return NULL;
    }
  public:
    compute(const T& f) : call(f) { } 
};


template<typename T = std::function<void()>>
void task_list(const std::vector<T>& jobs) {
  const size_t nt = jobs.size();
  tbb::task* dummy = new(tbb::task::allocate_root()) tbb::empty_task;
  dummy->set_ref_count(nt+1);

  tbb::task_list tl;
  for (size_t i = 0; i != nt; ++i) {
    tbb::task* childtask = new(dummy->allocate_child()) compute<T>(jobs[i]);
    tl.push_back(*childtask);
  }

  dummy->spawn_and_wait_for_all(tl);
  dummy->destroy(*dummy);
}

template<typename T = std::function<void()>>
void parallel_for(const std::vector<T>& jobs) {
  const size_t nt = jobs.size();
  static tbb::affinity_partitioner ap;
  //tbb::affinity_partitioner ap;
  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, nt),
    [&jobs](const tbb::blocked_range<size_t>& r) {
      const auto& rend = r.end(); 
      for(size_t i=r.begin(); i!=rend; ++i)
        jobs[i]();
    }, ap);
}

/*
template<typename T = std::function<void()>>
void limited_parallel_for(const std::vector<T>& jobs, const int max_num_threads) {
  tbb::task_arena limited_arena(max_num_threads);

  limited_arena.execute([&jobs] {
    const size_t nt = jobs.size();
    tbb::affinity_partitioner ap;
    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, nt),
      [&jobs](const tbb::blocked_range<size_t>& r) {
        const auto& rend = r.end(); 
        for(size_t i=r.begin(); i!=rend; ++i)
          jobs[i]();
      }, ap);
  });
}
*/

/*inline void parallel_for(const std::vector<std::function<void()>>& jobs) {
  tbb::parallel_for(
    tbb::blocked_range< const std::function<void()>* >(jobs.data(), jobs.data() + jobs.size()),
    [](const tbb::blocked_range< const std::function<void()>* >& r) {
      const auto& rend = r.end();
      for(auto it = r.begin(); it != rend; ++it)
        (*it)();
    });
}*/

template<typename T = std::function<void()>>
void parallel_for_each(const std::vector<T>& jobs) {
  tbb::parallel_for_each(jobs.begin(), jobs.end(), [](const std::function<void()>& call){ 
    call();
  });
}
#endif
