#include <gtest/gtest.h>

TEST(first, test) {
    EXPECT_STRNE("hello", "world");
    EXPECT_EQ(7*6, 42);
}
