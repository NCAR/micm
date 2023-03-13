#include <micm/configure/factory.hpp>

#include <gtest/gtest.h>
#include <string>
#include <iostream>

class NoArgs {
  public:
  std::string hello(){
    return "hello";
  }
};

TEST(Factory, DefaultConstructor){
  micm::Factory<NoArgs, std::string> factory;
}

TEST(Factory, CanCreateClassTakingNoArugments){
  micm::Factory<NoArgs, std::string> factory;

  factory.Register("NoArgs", [](){ return new NoArgs();});

  auto created = factory.CreateObject("NoArgs");
  EXPECT_EQ(created->hello(), "hello");
}

class StringArg {
  public:
    StringArg(std::string s)
      : s_(s) { }

    std::string hello(){ return s_;}
    std::string s_;
};
TEST(Factory, CanCreateClassTaking1StringArugment){
  using ObjectCreator = std::function<StringArg*(std::string)>;
  micm::Factory<StringArg, std::string, ObjectCreator> factory;

  factory.Register("StringArg", [](std::string arg){ return new StringArg(arg);});

  auto created = factory.CreateObject("StringArg", "hello");
  EXPECT_EQ(created->hello(), "hello");
}