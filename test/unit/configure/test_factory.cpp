#include <micm/configure/factory.hpp>

#include <gtest/gtest.h>
#include <string>
#include <iostream>

#ifdef USE_JSON
#include <nlohmann/json.hpp>
#endif

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

#ifdef USE_JSON
using json = nlohmann::json;
class JsonArg {
  public:
    std::string s_;
    int i_;
    double d_;

    JsonArg(const json& s)
      : s_(s["s"].get<std::string>())
      , i_(s["i"].get<int>())
      , d_(s["d"].get<double>())
      {
      }
};

TEST(Factory, CanCreateClassTakingJsonArgument){
  using ObjectCreator = std::function<JsonArg*(json)>;
  micm::Factory<JsonArg, std::string, ObjectCreator> factory;

  factory.Register("JsonArg", [](json arg){ return new JsonArg(arg);});

  json config = json::parse(R"(
    {
      "s": "hello",
      "i": 10,
      "d": 3.14
    }
  )");

  auto created = factory.CreateObject("JsonArg", config);
  EXPECT_EQ(created->s_, "hello");
  EXPECT_EQ(created->i_, 10);
  EXPECT_EQ(created->d_, 3.14);
}

// TEST(Factory, CanCreateMultipleClassesTakingJson){
//   using ObjectCreator = std::function<JsonArg*(json)>;
//   micm::Factory<JsonArg, std::string, ObjectCreator> factory;

//   factory.Register("JsonArg", [](json arg){ return new JsonArg(arg);});

//   json config = json::parse(R"(
//     {
//       "s": "hello",
//       "i": 10,
//       "d": 3.14
//     }
//   )");

//   auto created = factory.CreateObject("JsonArg", config);
//   EXPECT_EQ(created->s_, "hello");
//   EXPECT_EQ(created->i_, 10);
//   EXPECT_EQ(created->d_, 3.14);
// }
#endif