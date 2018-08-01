#ifndef VBL_H
#define VBL_H  // header guard

// to use new and old version of TimersAdvance and TimersAdvanceUntil simultaniously
#include <boost/optional.hpp>

//for timing output
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>      // std::setprecision  std::scientific

#include <boost/unordered_map.hpp>
#include <cmath>
#include <cstdio>

#include "vbl/CellsSystem.h"
#include "vbl/Environment.h"
#include "vbl/EnvironmentalSignals.h"

#include "vbl/geometry.h"

// namespace vbl
// {
//   vbl::CellsSystem *p_to_current_CellsSystem = NULL;
// }

//this pointer is necessary to keep track of the cells CellsSystem
// note: you need to initialize that when dealing with other software e.g. tumorcode!
#endif // header guard
