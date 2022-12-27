#include <micm/version.hpp>
#include <string>

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