#include <gtest/gtest.h>

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>

TEST(JitFunction, SimpleIntFunction)
{
  auto jit{micm::JitCompiler::create()};
  if (auto err = jit.takeError()) {
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