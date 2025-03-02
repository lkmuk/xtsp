#pragma once

#include <gtest/gtest.h>
#include <exception>
#include <functional>
#include <type_traits>
#include <string_view>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/bundled/format.h>

/// @brief More powerful than gtest's EXPECT_THROW.
/// Actually I prefer manually copying this so that EXPECT_* / ASSERT_* 
/// can pinpoint the specific source file, not this one.
template <typename ExceptionTy = std::invalid_argument>
inline void expect_throw_thisMsg(std::function<void()> faultyAction, const std::string& expectedMsg)
{

    static_assert(std::is_base_of<std::exception, ExceptionTy>::value);
    try
    {
        faultyAction();
        FAIL() << fmt::format("Should have thrown an exception: {}", expectedMsg);
    } 
    catch (const ExceptionTy& actualException)
    {
        EXPECT_EQ(actualException.what(), expectedMsg);
    }
}
    


