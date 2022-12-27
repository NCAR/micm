#include <micm/version.hpp>
#include <string>
#include <assert.h>

int main(){
  bool failed{false};

  try{
    std::string version = getmicmVersion();
  }
  catch(...)
  {
    failed = true;
  }

  assert(!failed);
}