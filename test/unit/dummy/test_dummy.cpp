#include <micm/dummy.hpp>
#include <assert.h>

int main(){
  bool failed{false};

  try{
    dummy();
  }
  catch(...)
  {
    failed = true;
  }

  assert(!failed);
}