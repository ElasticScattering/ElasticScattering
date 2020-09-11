#pragma once
#include <doctest.h>

#include "CPUComparisonTest.h"
#include "TestDetails.h"
#include "UtilKernelsTest.h"
#include "Test.h"

void test_main() {
    doctest::Context context;
    context.setOption("order-by", "file");
    context.setOption("no-breaks", true);
    int res = context.run();
}