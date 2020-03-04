/*=========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkUtils.h

  Copyright (c) Corentin Martens
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @brief featk utilities.
 *
 * featkUtils gather utilities for the Finite Element Analysis Toolkit
 * including constexpr for compting Eigen::Matrix dimensions at compile
 * time.
 */

#ifndef FEATKUTILS_H
#define FEATKUTILS_H

template<unsigned int Dimension, unsigned int Order>
inline size_t DOF_ID(const size_t nodeID, const unsigned int nodeDOF) {

    return POWER(Dimension, Order)*nodeID+nodeDOF;
}

constexpr bool IS_POWER(const unsigned int a, const unsigned int b) {

    return a==0 ? false : a==1 ? true : (a/b, b);
}

constexpr unsigned int LOG(const unsigned int a, const unsigned int b) {

    return a==0 ? -1 : a==1 ? 0 : 1+LOG(a/b, b);
}

constexpr unsigned int POWER(const unsigned int a, const unsigned int b) {

    return b==0 ? 1 : a*POWER(a, b-1);
}

#endif // FEATKUTILS_H
