/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Modified Factory class found in chapter 8 of
 * Alexandrescu, A., 2001. Modern C++ design: generic programming and design patterns applied. Addison-Wesley Longman
 * Publishing Co., Inc., USA.
 *
 */

#pragma once

#include <map>
#include <memory>

namespace micm
{

  template<class IdentifierType, class Object>
  class ThrowExceptionOnFactoryError
  {
   public:
    class Exception : public std::exception
    {
     public:
      const IdentifierType unknownId_;

     public:
      Exception(const IdentifierType& unknownId)
          : unknownId_(unknownId)
      {
      }

      virtual const char* what()
      {
        return "Unknown object type passed to Factory.";
      }
    };

   protected:
    static std::unique_ptr<Object> OnUnknownType(const IdentifierType& id)
    {
      throw Exception(id);
    }
  };

  template<
      class Object,
      typename IdentifierType,
      typename ObjectCreator = Object* (*)(),
      template<typename, class> class FactoryErrorPolicy = ThrowExceptionOnFactoryError>

  class Factory : public FactoryErrorPolicy<IdentifierType, Object>
  {
   public:
    bool Register(const IdentifierType& id, ObjectCreator creator)
    {
      return associations_.insert(typename AssocMap::value_type(id, creator)).second;
    }

    bool Unregister(const IdentifierType& id)
    {
      return associations_.erase(id) == 1;
    }

    template<typename... Args>
    std::unique_ptr<Object> CreateObject(const IdentifierType& id, Args&&... args)
    {
      auto i = associations_.find(id);
      if (i != associations_.end())
      {
        return std::unique_ptr<Object>((i->second)(std::forward<Args>(args)...));
      }
      else
      {
        return this->OnUnknownType(id);
      }
    }

   private:
    typedef std::map<IdentifierType, ObjectCreator> AssocMap;

    AssocMap associations_;
  };

}  // namespace micm