#include <gtest/gtest.h>

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>

TEST(JitFunction, SimpleIntFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo")
                               .arguments({ { "foo", micm::JitType::Int }, { "bar", micm::JitType::Int } })
                               .return_type(micm::JitType::Int);
  llvm::Value *ret_val = func.builder_->CreateNSWAdd(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  int (*func_ptr)(int, int) = (int (*)(int, int))(intptr_t)func_target.second;
  EXPECT_EQ(12, func_ptr(8, 4));
  EXPECT_EQ(-4, func_ptr(-8, 4));
  EXPECT_EQ(92, func_ptr(80, 12));
  func.exit_on_error_(func_target.first->remove());
}

TEST(JitFunction, SimpleFloatFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_float")
                               .arguments({ { "foo", micm::JitType::Float }, { "bar", micm::JitType::Float } })
                               .return_type(micm::JitType::Float);
  llvm::Value *ret_val = func.builder_->CreateFAdd(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  float (*func_ptr)(float, float) = (float (*)(float, float))(intptr_t)func_target.second;
  EXPECT_EQ(8.32f + 4.23f, func_ptr(8.32f, 4.23f));
  EXPECT_EQ(-8.93f + 4.01f, func_ptr(-8.93f, 4.01f));
  EXPECT_EQ(80.12f + 12.42f, func_ptr(80.12f, 12.42f));
  func.exit_on_error_(func_target.first->remove());
}

TEST(JitFunction, SimpleDoubleFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_double")
                               .arguments({ { "foo", micm::JitType::Double }, { "bar", micm::JitType::Double } })
                               .return_type(micm::JitType::Double);
  llvm::Value *ret_val = func.builder_->CreateFAdd(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  double (*func_ptr)(double, double) = (double (*)(double, double))(intptr_t)func_target.second;
  EXPECT_EQ(8.32 + 4.23, func_ptr(8.32, 4.23));
  EXPECT_EQ(-8.93 + 4.01, func_ptr(-8.93, 4.01));
  EXPECT_EQ(80.12 + 12.42, func_ptr(80.12, 12.42));
  func.exit_on_error_(func_target.first->remove());
}

TEST(JitFunction, SimpleBoolFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_bool")
                               .arguments({ { "foo", micm::JitType::Bool }, { "bar", micm::JitType::Bool } })
                               .return_type(micm::JitType::Bool);
  llvm::Value *ret_val = func.builder_->CreateOr(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  bool (*func_ptr)(bool, bool) = (bool (*)(bool, bool))(intptr_t)func_target.second;
  EXPECT_EQ(true, func_ptr(true, false));
  EXPECT_EQ(true, func_ptr(true, true));
  EXPECT_EQ(false, func_ptr(false, false));
  func.exit_on_error_(func_target.first->remove());
}