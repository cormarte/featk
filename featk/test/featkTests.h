#ifndef FEATKTESTS_H
#define FEATKTESTS_H

#include <featk/core/featkGlobal.h>

#define EPS 1.0E-4

FEATK_EXPORT bool featkHex8StiffnessMatrixTest();
FEATK_EXPORT void featkRunAllTests();
FEATK_EXPORT bool featkTet4StiffnessMatrixTest();
FEATK_EXPORT bool featkTet4LinearElasticitySolverTest();

#endif // FEATKTESTS_H
