#include <gtest/gtest.h>
#include "../static.h"
#ifdef HAVE_TBB
  #include "../thread/tbb_interface.h"
#endif

int main(int argc, char **argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  static_variables();

#ifdef HAVE_TBB
  tbb::task_scheduler_init init(tbb::task_scheduler_init::automatic);
  //tbb::task_scheduler_init init(1);
#endif
 
  result = RUN_ALL_TESTS();
  return result;
}
