//
// Created by sbian on 2019/9/3.
//

#ifndef BIM_STDAFX_H
#define BIM_STDAFX_H

#include "targetver.h"

#include <stdio.h>
#ifdef _WIN32
#include <tchar.h>
#endif

// TODO: reference additional headers your program requires here
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <ctime>
#include <memory>
#include <random>
#include <functional>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <tuple>
#include <deque>

#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
#define DSFMT_MEXP 19937
#endif
#endif
#include "SFMT/dSFMT/dSFMT.h"
#include "timer.h"
#include "commonStruct.h"
#include "commonFunc.h"
#include "serialize.h"
#include "resultInfo.h"
#include "IOcontroller.h"
#include "argument.h"
#include "graphBase.h"
#include "hyperGraph.h"

#include "celf.h"
#include "celfpp.h"
#include "alg.h"

#endif //BIM_STDAFX_H
