#include <gtest/gtest.h>

#include <memory>
#include <interface.hpp>

TEST(FortranInterface, CanCallCPPFunction){
  micm::get_solver("filepath.txt");
}

TEST(FortranInterface, CanCallSolver){
  micm::FuncPtr f = micm::get_solver("filepath.txt");
  std::unique_ptr<double[]> arg1 = std::make_unique<double[]>(5);
  std::unique_ptr<double[]> arg2 = std::make_unique<double[]>(5);
  std::unique_ptr<double[]> arg3 = std::make_unique<double[]>(5);
  f(arg1.get(), arg2.get(), arg3.get());
}