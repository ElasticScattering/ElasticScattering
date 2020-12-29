#pragma once
#include <doctest.h>

#include "TestMacros.h"

#include "GridGenerationTests.h"

#include "AngleTests.h"
#include "OrbitDetailsTest.h"

#include "SimpsonWeightsTest.h"

#include "IndexTests.h"
#include "GridMovementTests.h"

#include "LifetimeTests.h"
#include "IntegrationTest.h"

//#include "BasicOpenCLTest.h"
//#include "UtilKernelsTest.h"

void test_main() {
    doctest::Context context;
    context.setOption("order-by", "file");
    context.setOption("no-breaks", true);
    //context.setOption("success", true);
    int res = context.run();
}
