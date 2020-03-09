/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkTet4Element.h

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
 * @class featkTet4Element
 *
 * @brief Specialization of a 4-node tetrahedral featkElement.
 *
 * featkTet4Element is a specialization of featkElement representing an
 * isoparametric tetrahedron in Dimension 3 with 4 Nodes, 4 Boundaries
 * parametrized in NaturalDimension 3. featkTet4Element provides
 * implementations for the specialization of member variables
 * featkElement::integrationRule and featkElement::nodesNaturalCoordinates
 * and member functions featkElement::featkElement(),
 * featkElement::getShapeFunctionNaturalDerivativeValues() and
 * featkElement::getShapeFunctionValues().
 *
 * featkTet4Element is parametrized as follows:
 *
 *
 *     (3)
 *      |\ \(2)
 *      | \ //\
 *      |  /\//\         ζ
 *      | ///\//\        | η
 *      |//////\/\       |/
 *     (0)-------(1)     o----ξ
 *
 *
 * |  Node  |  ξ  |  η  |  ζ  |
 * | :----: | :-: | :-: | :-: |
 * |  (0)   |  0  |  0  |  0  |
 * |  (1)   |  1  |  0  |  0  |
 * |  (2)   |  0  |  1  |  0  |
 * |  (3)   |  0  |  0  |  1  |
 *
 *
 * \f[
 *
 * \begin{align}
 * N_0(\xi, \eta, \zeta) &= 1−\xi-\eta-\zeta\\
 * N_1(\xi, \eta, \zeta) &= \xi\\
 * N_2(\xi, \eta, \zeta) &= \eta\\
 * N_3(\xi, \eta, \zeta) &= \zeta\\
 * \end{align}
 *
 * \f]
 *
 *
 * See J. E. Flaherty, Finite Element Analysis: Chapter 4 - Finite Element
 * Approximation, p.24, Figure 4.5.2.
 *
 */

#ifndef FEATKTET4ELEMENT_H
#define FEATKTET4ELEMENT_H

#include <featk/geometry/featkElement.h>

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension> class featkElement;

using featkTet4Element = featkElement<3, 4, 4, 3>;

#endif // FEATKTET4ELEMENT_H
