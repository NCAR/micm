#include <micm/dummy.hpp>
#include <assert.h>

int main(){
  bool failed{false};

  try{
    double d = dummy(1.1, 2.2);
    decl();
  }
  catch(...)
  {
    failed = true;
  }

  assert(!failed);
}