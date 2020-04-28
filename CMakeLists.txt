cmake_minimum_required(VERSION 3.10)
project(BIM)

set(CMAKE_CXX_STANDARD 14)

add_executable(BIM main.cpp targetver.h stdafx.h argument.h commonStruct.h graphBase.h commonFunc.h IOcontroller.h serialize.h resultInfo.h SFMT/dSFMT/dSFMT-common.h SFMT/dSFMT/dSFMT-params.h SFMT/dSFMT/dSFMT-params521.h SFMT/dSFMT/dSFMT-params1279.h SFMT/dSFMT/dSFMT-params2203.h SFMT/dSFMT/dSFMT-params4253.h SFMT/dSFMT/dSFMT-params11213.h SFMT/dSFMT/dSFMT-params19937.h SFMT/dSFMT/dSFMT-params44497.h SFMT/dSFMT/dSFMT-params86243.h SFMT/dSFMT/dSFMT-params132049.h SFMT/dSFMT/dSFMT-params216091.h SFMT/dSFMT/dSFMT.h SFMT/dSFMT/test.c timer.h alg.h hyperGraph.h celf.h celfpp.h)