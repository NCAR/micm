#include <iostream>
#include <string>
#include <vector>

class Test {
  public:
    void function(auto X, auto Y);
};

void Test::function(auto X, auto Y)
{
  for(const auto& x : X)
  {
    std::cout << x << " ";
  }
  std::cout << std::endl;
  for(const auto& y : Y)
  {
    std::cout << y << " ";
  }
  std::cout << std::endl;
}

int main()
{
  Test t;
  std::vector<double> vd = {1.0, 2.0, 3.0};
  std::vector<int> vi = {1, 2, 3};
  std::vector<float> vf = {1.0, 2.0, 3.0};
  std::vector<std::string> vs = {"a", "b", "c"};

  t.function(vd, vd);
  t.function(vi, vi);
  t.function(vf, vf);
  t.function(vs, vs);
  t.function(vd, vi);
  t.function(vi, vf);
  t.function(vf, vs);
  t.function(vs, vd);
  return 0;
}