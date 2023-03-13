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

  template<class Object, typename IdentifierType, typename ObjectCreator, template<typename, class> class FactoryErrorPolicy>

  class Factory : public FactoryErrorPolicy<IdentifierType, Object>
  {
    public:
      bool Register(const IdentifierType& id, ObjectCreator creator) {
        return associations_.insert(
          AssocMap::value_type(id, creator)
        ).second;
      }

      bool Unregister(const IdentifierType& id)
      {
        return associations_.erase(id) == 1;
      }

      std::unique_ptr<Object> CreateObject(const IdentifierType& id)
      {
        auto i = associations_.find(id);
        if (i != associations_.end())
        {
          return (i->second)();
        }
        else {
          // TODO: what kind of error should we throw?
          // In general, anything in the generics/configure directory (which probably should be combined?)
          // I think should throw errors because an incorrectly configured system cannot run
        }
      }

    private:
      typedef std::map<IdentifierType, ObjectCreator> AssocMap;

      AssocMap associations_;
  };

}  // namespace micm