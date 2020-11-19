#pragma once
#include <doctest.h>

#include "OrbitDetailsTest.h"
#include "ImpurityGridTests.h"
#include "SimpsonWeightsTest.h"

//#include "BasicOpenCLTest.h"
//#include "UtilKernelsTest.h"

void test_main() {
    doctest::Context context;
    context.setOption("order-by", "file");
    context.setOption("no-breaks", true);
    int res = context.run();
}
