// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <memory>  // add this include

template<class T>
class ScalarView
{
    std::shared_ptr<T> data_;

public:
    using value_type = T;

    ScalarView(T init = T{}) : data_(std::make_shared<T>(init)) {}

    /// Returns T& even from a const ScalarView — the shared_ptr can't be
    /// reseated but the pointed-to value is still mutable. This is what
    /// lets [=]-captured copies accumulate into shared storage.
    MICM_DEVICE_FUNCTION T& host_value() const { return *data_; }

    operator T() const { return *data_; }
    T& operator=(const T& val) { *data_ = val; return *data_; }
    ScalarView& operator=(const ScalarView<T>& val) { *data_ = *(val.data_); return *this; }
    T* data() const // NOLINT(readability-identifier-naming)
    {
        return data_.get();
    }

    void CopyToHost() const {}
    void CopyToDevice() const {}
};