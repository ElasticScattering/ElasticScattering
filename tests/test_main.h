#pragma once
#include <doctest.h>

#include "TestMacros.h"

#include "es/AngleTests.h"
#include "es/OrbitDetailsTest.h"
#include "es/IndexTests.h"
#include "es/GridMovementTests.h"
#include "es/LifetimeTests.h"

#include "sim/GridGenerationTests.h"
#include "sim/SimpsonWeightsTest.h"
#include "sim/SigmaIntegrationTest.h"

//#include "sim/cl/BasicOpenCLTest.h"
//#include "UtilKernelsTest.h"

void test_main() {
    doctest::Context context;
    context.setOption("order-by", "file");
    context.setOption("no-breaks", true);
    //context.setOption("success", true);
    int res = context.run();
}
