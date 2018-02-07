#include <memory>
#include <fstream>
#include <iostream>
#include "tag_generator.h"
#include "shared_tensor.h"

using namespace std;

TAG_Generator* tag__;
std::unique_ptr<TAG_Generator> tag;

int main() {
  tag = unique_ptr<TAG_Generator>(new TAG_Generator());
  tag__ = tag.get();
  
  //for (size_t i = 0; i != 10; ++i) {
  //  string filename = tag__->get(); 
  //  std::ofstream o(filename);
  //}
 
  
  Tensor tensor({3u,1u,2u});
  tensor.randomize();
  cout << "Tensor:\n" << tensor << endl;
 
  Shared_Tensor obj(make_shared<Tensor>(tensor), true);  
  
  cout << "COUNT = " << obj.use_count() << endl;
  
  obj.load();
  cout << "COUNT = " << obj.use_count() << endl;
  obj.access().scale(100);
  obj.dump();
  cout << "COUNT = " << obj.use_count() << endl;

  //const Shared_Tensor obj2(obj);
  //cout << "COUNT = " << obj.use_count() << " :: " << obj2.use_count() << endl;
  //obj2.load();
  //cout << "Tensor:\n" << obj2.access() << endl;
  //cout << "COUNT = " << obj.use_count() << " :: " << obj2.use_count() << endl;
  //obj2.dump();
  //cout << "COUNT = " << obj.use_count() << " :: " << obj2.use_count() << endl;

  Tensor tensor2({3u,1u,2u});
  auto&& ref = obj; 
  ref.load();
  cout << "COUNT = " << obj.use_count() << endl;
  ref.access() = tensor2;
  ref.dump();
  cout << "COUNT = " << obj.use_count() << endl;

  obj.load();
  cout << "Tensor:\n" << obj.access() << endl;
  obj.dump();

  tag__->reset();
  return 0;
}
