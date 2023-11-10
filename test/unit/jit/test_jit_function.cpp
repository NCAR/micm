#include <gtest/gtest.h>

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/jit_function.hpp>

// This test creates a function that adds two integers and returns the sum
TEST(JitFunction, SimpleInt32Function)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_int32")
                               .arguments({ { "foo", micm::JitType::Int32 }, { "bar", micm::JitType::Int32 } })
                               .return_type(micm::JitType::Int32);
  llvm::Value *ret_val = func.builder_->CreateNSWAdd(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);

  auto func_target = func.Generate();
  int32_t (*func_ptr)(int32_t, int32_t) = (int32_t(*)(int32_t, int32_t))(intptr_t)func_target.second;
  EXPECT_EQ(12, func_ptr(8, 4));
  EXPECT_EQ(-4, func_ptr(-8, 4));
  EXPECT_EQ(92, func_ptr(80, 12));
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that adds two integers and returns the sum
TEST(JitFunction, SimpleInt64Function)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_int64")
                               .arguments({ { "foo", micm::JitType::Int64 }, { "bar", micm::JitType::Int64 } })
                               .return_type(micm::JitType::Int64);
  llvm::Value *ret_val = func.builder_->CreateNSWAdd(func.arguments_[0].ptr_, func.arguments_[1].ptr_, "add args");
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  int64_t (*func_ptr)(int64_t, int64_t) = (int64_t(*)(int64_t, int64_t))(intptr_t)func_target.second;
  EXPECT_EQ(12l, func_ptr(8l, 4l));
  EXPECT_EQ(-4l, func_ptr(-8l, 4l));
  EXPECT_EQ(92l, func_ptr(80l, 12l));
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that adds two floats and returns the sum
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

// This test creates a function that adds two doubles and returns the sum
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

// This test creates a function that returns an OR of two booleans
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

// This test creates a function that adds the third elements in two integer arrays, sets the
// second element of the second array as the sum, and also returns the sum
TEST(JitFunction, SimpleInt32PtrFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_int32_ptr")
                               .arguments({ { "foo", micm::JitType::Int32Ptr }, { "bar", micm::JitType::Int32Ptr } })
                               .return_type(micm::JitType::Int32);
  llvm::Value *index_list[1];
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 2));
  llvm::Value *foo_val = func.GetArrayElement(func.arguments_[0], index_list, micm::JitType::Int32);
  llvm::Value *bar_val = func.GetArrayElement(func.arguments_[1], index_list, micm::JitType::Int32);
  llvm::Value *sum = func.builder_->CreateNSWAdd(foo_val, bar_val, "sum foo bar");
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1));
  func.SetArrayElement(func.arguments_[1], index_list, micm::JitType::Int32, sum);
  func.builder_->CreateRet(sum);
  auto func_target = func.Generate();
  int32_t (*func_ptr)(int32_t *, int32_t *) = (int32_t(*)(int32_t *, int32_t *))(intptr_t)func_target.second;
  int32_t a[] = { 9, 4, 33 };
  int32_t b[] = { 4, 21, 2, 42 };
  EXPECT_EQ(35, func_ptr(a, b));
  EXPECT_EQ(9, a[0]);
  EXPECT_EQ(4, a[1]);
  EXPECT_EQ(33, a[2]);
  EXPECT_EQ(4, b[0]);
  EXPECT_EQ(35, b[1]);
  EXPECT_EQ(2, b[2]);
  EXPECT_EQ(42, b[3]);
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that adds the third elements in two integer arrays, sets the
// second element of the second array as the sum, and also returns the sum
TEST(JitFunction, SimpleInt64PtrFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_int64_ptr")
                               .arguments({ { "foo", micm::JitType::Int64Ptr }, { "bar", micm::JitType::Int64Ptr } })
                               .return_type(micm::JitType::Int64);
  llvm::Value *index_list[1];
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 2));
  llvm::Value *foo_val = func.GetArrayElement(func.arguments_[0], index_list, micm::JitType::Int64);
  llvm::Value *bar_val = func.GetArrayElement(func.arguments_[1], index_list, micm::JitType::Int64);
  llvm::Value *sum = func.builder_->CreateNSWAdd(foo_val, bar_val, "sum foo bar");
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1));
  func.SetArrayElement(func.arguments_[1], index_list, micm::JitType::Int64, sum);
  func.builder_->CreateRet(sum);
  auto func_target = func.Generate();
  int64_t (*func_ptr)(int64_t *, int64_t *) = (int64_t(*)(int64_t *, int64_t *))(intptr_t)func_target.second;
  int64_t a[] = { 9l, 4l, 33l };
  int64_t b[] = { 4l, 21l, 2l, 42l };
  EXPECT_EQ(35l, func_ptr(a, b));
  EXPECT_EQ(9l, a[0]);
  EXPECT_EQ(4l, a[1]);
  EXPECT_EQ(33l, a[2]);
  EXPECT_EQ(4l, b[0]);
  EXPECT_EQ(35l, b[1]);
  EXPECT_EQ(2l, b[2]);
  EXPECT_EQ(42l, b[3]);
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that adds the third elements in two float arrays, sets the
// second element of the second array as the sum, and also returns the sum
TEST(JitFunction, SimpleFloatPtrFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_float_ptr")
                               .arguments({ { "foo", micm::JitType::FloatPtr }, { "bar", micm::JitType::FloatPtr } })
                               .return_type(micm::JitType::Float);
  llvm::Value *index_list[1];
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 2));
  llvm::Value *foo_val = func.GetArrayElement(func.arguments_[0], index_list, micm::JitType::Float);
  llvm::Value *bar_val = func.GetArrayElement(func.arguments_[1], index_list, micm::JitType::Float);
  llvm::Value *sum = func.builder_->CreateFAdd(foo_val, bar_val, "sum foo bar");
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1));
  func.SetArrayElement(func.arguments_[1], index_list, micm::JitType::Float, sum);
  func.builder_->CreateRet(sum);
  auto func_target = func.Generate();
  float (*func_ptr)(float *, float *) = (float (*)(float *, float *))(intptr_t)func_target.second;
  float a[] = { 9.3f, 4.4f, 33.3f };
  float b[] = { 4.53f, 21.02f, 2.0f, 42.23f };
  EXPECT_EQ(35.3f, func_ptr(a, b));
  EXPECT_EQ(9.3f, a[0]);
  EXPECT_EQ(4.4f, a[1]);
  EXPECT_EQ(33.3f, a[2]);
  EXPECT_EQ(4.53f, b[0]);
  EXPECT_EQ(35.3f, b[1]);
  EXPECT_EQ(2.0f, b[2]);
  EXPECT_EQ(42.23f, b[3]);
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that adds the third elements in two double arrays, sets the
// second element of the second array as the sum, and also returns the sum
TEST(JitFunction, SimpleDoublePtrFunction)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo_double_ptr")
                               .arguments({ { "foo", micm::JitType::DoublePtr }, { "bar", micm::JitType::DoublePtr } })
                               .return_type(micm::JitType::Double);
  llvm::Value *index_list[1];
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 2));
  llvm::Value *foo_val = func.GetArrayElement(func.arguments_[0], index_list, micm::JitType::Double);
  llvm::Value *bar_val = func.GetArrayElement(func.arguments_[1], index_list, micm::JitType::Double);
  llvm::Value *sum = func.builder_->CreateFAdd(foo_val, bar_val, "sum foo bar");
  index_list[0] = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1));
  func.SetArrayElement(func.arguments_[1], index_list, micm::JitType::Double, sum);
  func.builder_->CreateRet(sum);
  auto func_target = func.Generate();
  double (*func_ptr)(double *, double *) = (double (*)(double *, double *))(intptr_t)func_target.second;
  double a[] = { 9.3, 4.4, 33.3 };
  double b[] = { 4.53, 21.02, 2.0, 42.23 };
  EXPECT_EQ(35.3, func_ptr(a, b));
  EXPECT_EQ(9.3, a[0]);
  EXPECT_EQ(4.4, a[1]);
  EXPECT_EQ(33.3, a[2]);
  EXPECT_EQ(4.53, b[0]);
  EXPECT_EQ(35.3, b[1]);
  EXPECT_EQ(2.0, b[2]);
  EXPECT_EQ(42.23, b[3]);
  func.exit_on_error_(func_target.first->remove());
}

// This test creates a function that includes a loop that adds 1 to a variable that starts at 0
// and iterates 10 times and returns the summed value
TEST(JitFunction, SimpleLoop)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func =
      micm::JitFunction::create(jit.get()).name("foo_loop").arguments({}).return_type(micm::JitType::Int32);
  auto loop = func.StartLoop("foo loop", 0, 10);
  llvm::PHINode *ret_val = func.builder_->CreatePHI(func.GetType(micm::JitType::Int32), 2, "ret val");
  ret_val->addIncoming(llvm::ConstantInt::get(*(func.context_), llvm::APInt(32, 1)), func.entry_block_);
  llvm::Value *incr = llvm::ConstantInt::get(*(func.context_), llvm::APInt(32, 1));
  llvm::Value *next_val = func.builder_->CreateNSWAdd(ret_val, incr, "add incr");
  func.EndLoop(loop);
  ret_val->addIncoming(next_val, loop.block_);
  func.builder_->CreateRet(ret_val);
  auto func_target = func.Generate();
  int32_t (*func_ptr)() = (int32_t(*)())(intptr_t)func_target.second;
  EXPECT_EQ(10, func_ptr());
  func.exit_on_error_(func_target.first->remove());
}

// This test creates two functions, foo and bar, that return the sum of their single integer
// arguments and 10 and 100, respectively
TEST(JitFunction, MultipleFunctions)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction foo_func = micm::JitFunction::create(jit.get())
                                   .name("foo")
                                   .arguments({ { "arg", micm::JitType::Int32 } })
                                   .return_type(micm::JitType::Int32);
  llvm::Value *foo_ret_val = foo_func.builder_->CreateNSWAdd(
      foo_func.arguments_[0].ptr_, llvm::ConstantInt::get(*(foo_func.context_), llvm::APInt(32, 10)), "add args");
  foo_func.builder_->CreateRet(foo_ret_val);
  auto foo_target = foo_func.Generate();
  int32_t (*foo_func_ptr)(int) = (int32_t(*)(int))(intptr_t)foo_target.second;
  micm::JitFunction bar_func = micm::JitFunction::create(jit.get())
                                   .name("bar")
                                   .arguments({ { "arg", micm::JitType::Int32 } })
                                   .return_type(micm::JitType::Int32);
  llvm::Value *bar_ret_val = bar_func.builder_->CreateNSWAdd(
      bar_func.arguments_[0].ptr_, llvm::ConstantInt::get(*(bar_func.context_), llvm::APInt(32, 100)), "add args");
  bar_func.builder_->CreateRet(bar_ret_val);
  auto bar_target = bar_func.Generate();
  int32_t (*bar_func_ptr)(int) = (int32_t(*)(int))(intptr_t)bar_target.second;
  EXPECT_EQ(32, foo_func_ptr(22));
  EXPECT_EQ(102, bar_func_ptr(2));
  EXPECT_EQ(254, bar_func_ptr(foo_func_ptr(144)));
  foo_func.exit_on_error_(foo_target.first->remove());
  bar_func.exit_on_error_(bar_target.first->remove());
}

// This test creates a local array of ints, populates it with a specified value,
// and returns the sum
TEST(JitFunction, LocalArray)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func = micm::JitFunction::create(jit.get())
                               .name("foo")
                               .arguments({ { "arg", micm::JitType::Int64 } })
                               .return_type(micm::JitType::Int64);

  auto int_type = func.GetType(micm::JitType::Int64);
  llvm::Value *zero = llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 0));
  llvm::Type *foo_array_type = llvm::ArrayType::get(int_type, 10);
  llvm::AllocaInst *foo_array =
      func.builder_->CreateAlloca(foo_array_type, llvm::ConstantInt::get(*(func.context_), llvm::APInt(64, 1)), "foo_array");

  // loop to set array elements
  auto loop = func.StartLoop("set_loop", 0, 10);
  llvm::Value *index_list[2];
  index_list[0] = zero;
  index_list[1] = loop.index_;
  llvm::Value *set_elem = func.builder_->CreateInBoundsGEP(foo_array_type, foo_array, index_list, "set_elem_ptr");
  func.builder_->CreateStore(func.arguments_[0].ptr_, set_elem);
  func.EndLoop(loop);

  // loop to sum array elements
  index_list[1] = zero;
  llvm::Value *get_elem = func.builder_->CreateInBoundsGEP(foo_array_type, foo_array, index_list, "get_first_elem_ptr");
  llvm::Value *first_elem = func.builder_->CreateLoad(int_type, get_elem, "load_first_elem");
  loop = func.StartLoop("sum_loop", 0, 10, 2);
  llvm::PHINode *ret_val = func.builder_->CreatePHI(int_type, 2, "ret_val");
  ret_val->addIncoming(first_elem, loop.prior_block_);
  index_list[1] = loop.index_;
  get_elem = func.builder_->CreateInBoundsGEP(foo_array_type, foo_array, index_list, "get_curr_elem_ptr");
  llvm::Value *curr_elem = func.builder_->CreateLoad(int_type, get_elem, "load_curr_elem");
  llvm::Value *next_val = func.builder_->CreateNSWAdd(ret_val, curr_elem, "add_curr_elem");
  func.EndLoop(loop);
  ret_val->addIncoming(next_val, loop.block_);

  func.builder_->CreateRet(ret_val);

  auto foo_target = func.Generate();
  int64_t (*func_ptr)(int) = (int64_t(*)(int))(intptr_t)foo_target.second;
  EXPECT_EQ(20, func_ptr(4));
  EXPECT_EQ(5, func_ptr(1));
  func.exit_on_error_(foo_target.first->remove());
}

// This test creates several functions with the same name that
// add a unique number to the function argument and return the sum
TEST(JitFunction, SameNameFunctions)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error] ");
    EXPECT_TRUE(false);
  }
  micm::JitFunction func1 = micm::JitFunction::create(jit.get())
                                .name("foobar")
                                .arguments({ { "foo", micm::JitType::Int32 } })
                                .return_type(micm::JitType::Int32);
  llvm::Value *const_val = llvm::ConstantInt::get(*(func1.context_), llvm::APInt(64, 2));
  llvm::Value *ret_val = func1.builder_->CreateNSWAdd(func1.arguments_[0].ptr_, const_val, "add args");
  func1.builder_->CreateRet(ret_val);

  auto func1_target = func1.Generate();
  int32_t (*func1_ptr)(int32_t) = (int32_t(*)(int32_t))(intptr_t)func1_target.second;

  micm::JitFunction func2 = micm::JitFunction::create(jit.get())
                                .name("foobar")
                                .arguments({ { "foo", micm::JitType::Int32 } })
                                .return_type(micm::JitType::Int32);
  const_val = llvm::ConstantInt::get(*(func2.context_), llvm::APInt(64, 12));
  ret_val = func2.builder_->CreateNSWAdd(func2.arguments_[0].ptr_, const_val, "add args");
  func2.builder_->CreateRet(ret_val);

  auto func2_target = func2.Generate();
  int32_t (*func2_ptr)(int32_t) = (int32_t(*)(int32_t))(intptr_t)func2_target.second;

  micm::JitFunction func3 = micm::JitFunction::create(jit.get())
                                .name("foobar")
                                .arguments({ { "foo", micm::JitType::Int32 } })
                                .return_type(micm::JitType::Int32);
  const_val = llvm::ConstantInt::get(*(func3.context_), llvm::APInt(64, 24));
  ret_val = func3.builder_->CreateNSWAdd(func3.arguments_[0].ptr_, const_val, "add args");
  func3.builder_->CreateRet(ret_val);

  auto func3_target = func3.Generate();
  int32_t (*func3_ptr)(int32_t) = (int32_t(*)(int32_t))(intptr_t)func3_target.second;

  EXPECT_EQ(10, func1_ptr(8));
  EXPECT_EQ(-6, func1_ptr(-8));
  EXPECT_EQ(82, func1_ptr(80));

  EXPECT_EQ(20, func2_ptr(8));
  EXPECT_EQ(4, func2_ptr(-8));
  EXPECT_EQ(92, func2_ptr(80));

  EXPECT_EQ(32, func3_ptr(8));
  EXPECT_EQ(16, func3_ptr(-8));
  EXPECT_EQ(104, func3_ptr(80));

  func1.exit_on_error_(func1_target.first->remove());
  func2.exit_on_error_(func2_target.first->remove());
  func3.exit_on_error_(func3_target.first->remove());
}