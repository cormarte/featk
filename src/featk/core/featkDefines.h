/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkDefines.h

  Copyright (c) Corentin Martens
  All rights reserved.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
     NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR
     ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR
     OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING
     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
     OTHER DEALINGS IN THE SOFTWARE.

==========================================================================*/

/**
 *
 * @brief Type definitions for the Finite Element Analysis Toolkit.
 *
 * featkDefines gather element and matrix type definitions for the Finite
 * Element Analysis Toolkit.
 *
 * featkElementType enumerated type is used by featkElementInterface to
 * identifiy its instantiated concrete featkElement type at run time.
 *
 * As the dimensions of the various matrices involved in finite element
 * problems are known at compile time given the cartesian dimension of the
 * problem, the number of nodes and natural dimension of the element types
 * involved and/or the order of the problem variable, featk can benefit from
 * Eigen's compile time optimization. featkDefines defines templated aliases
 * for these various matrices or better user readibility.
 *
 */

#ifndef FEATKDEFINES_H
#define FEATKDEFINES_H

#include <featk/core/featkUtils.h>

#include <Eigen/Dense>

using namespace Eigen;

enum featkElementType : unsigned char {FEATK_TET4, FEATK_HEX8};

template<unsigned int Dimension, unsigned int Order> using AttributeValueType = Matrix<double, POWER(Dimension, Order/2+Order%2), POWER(Dimension, Order/2)>;

template<unsigned int Dimension, unsigned int Nodes, unsigned int Order> using BMatrixType = Matrix<double, POWER(Dimension, Order+1), Nodes*POWER(Dimension, Order)>;
template<unsigned int Dimension, unsigned int Order>                     using CMatrixType = AttributeValueType<Dimension, Order>;
template<unsigned int Dimension, unsigned int Order>                     using DMatrixType = Matrix<double, POWER(Dimension, Order+1), 1>;
template<unsigned int Dimension, unsigned int Nodes, unsigned int Order> using KMatrixType = Matrix<double, Nodes*POWER(Dimension, Order), Nodes*POWER(Dimension, Order)>;
template<unsigned int Dimension, unsigned int Nodes, unsigned int Order> using NMatrixType = Matrix<double, POWER(Dimension, Order), Nodes*POWER(Dimension, Order)>;
template<unsigned int Dimension, unsigned int Nodes, unsigned int Order> using QMatrixType = Matrix<double, Nodes*POWER(Dimension, Order), 1>;

template<unsigned int Dimension, unsigned int NaturalDimension> using JacobianMatrixType = Matrix<double, NaturalDimension, Dimension>;
template<unsigned int NaturalDimension>                         using NaturalCoordinatesMatrixType = Matrix<double, 1, NaturalDimension>;
template<unsigned int Dimension, unsigned int Nodes>            using NodesCartesianCoordinatesMatrixType = Matrix<double, Nodes, Dimension>;
template<unsigned int Nodes, unsigned int NaturalDimension>     using NodesNaturalCoordinatesMatrixType = Matrix<double, Nodes, NaturalDimension>;
template<unsigned int Dimension, unsigned int Nodes>            using ShapeFunctionCartesianDerivativeValuesMatrixType = Matrix<double, Dimension, Nodes>;
template<unsigned int Nodes, unsigned int NaturalDimension>     using ShapeFunctionNaturalDerivativeValuesMatrixType = Matrix<double, NaturalDimension, Nodes>;
template<unsigned int Nodes>                                    using ShapeFunctionValuesMatrixType = Matrix<double, 1, Nodes>;

#endif // FEATKDEFINES_H
