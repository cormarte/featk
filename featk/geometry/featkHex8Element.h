/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkHex8Element.h

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
 * @class featkHex8Element
 *
 * @brief Specialization of a 8-node hexahedral featkElement.
 *
 * featkHex8Element is a specialization of featkElement representing an
 * isoparametric tetrahedron in Dimension 3 with 8 Nodes, 6 Boundaries
 * parametrized in NaturalDimension 3. featkHex8Element provides
 * implementations for the specialization of member variables
 * featkElement::integrationRule and featkElement::nodesNaturalCoordinates
 * and member functions featkElement::featkElement(),
 * featkElement::getShapeFunctionNaturalDerivativeValues() and
 * featkElement::getShapeFunctionValues().
 *
 * featkHex8Element is parametrization used as follows:
 *
 *
 *         (7)-------(6)
 *         /|        /|
 *        / |       / |      ζ
 *       /  |      /  |      | η
 *     (4)-(3)---(5)-(2)     |/
 *      |  ///////|///       o----ξ
 *      | ////////|//
 *      |/////////|/
 *     (0)-------(1)
 *
 *
 * |  Node  |  ξ  |  η  |  ζ  |
 * | :----: | :-: | :-: | :-: |
 * |  (0)   | −1  | −1  | −1  |
 * |  (1)   | +1  | −1  | −1  |
 * |  (2)   | +1  | +1  | −1  |
 * |  (3)   | −1  | +1  | −1  |
 * |  (4)   | −1  | −1  | +1  |
 * |  (5)   | +1  | −1  | +1  |
 * |  (6)   | +1  | +1  | +1  |
 * |  (7)   | −1  | +1  | +1  |
 *
 *
 *  N0(ξ, η, ζ) = (1/8)*(1−ξ)*(1−η)*(1−ζ)\n
 *  N1(ξ, η, ζ) = (1/8)*(1+ξ)*(1−η)*(1−ζ)\n
 *  N2(ξ, η, ζ) = (1/8)*(1+ξ)*(1+η)*(1−ζ)\n
 *  N3(ξ, η, ζ) = (1/8)*(1−ξ)*(1+η)*(1−ζ)\n
 *  N4(ξ, η, ζ) = (1/8)*(1−ξ)*(1−η)*(1+ζ)\n
 *  N5(ξ, η, ζ) = (1/8)*(1+ξ)*(1−η)*(1+ζ)\n
 *  N6(ξ, η, ζ) = (1/8)*(1+ξ)*(1+η)*(1+ζ)\n
 *  N7(ξ, η, ζ) = (1/8)*(1−ξ)*(1+η)*(1+ζ)\n
 *
 *
 * See C. Felippa, Advanced Finite Element Methods: Chapter 11 - The 8-Node
 * Hexahedron, p. 4, Tables 11.1-2, 2017.
 *
 */

#ifndef FEATKHEX8ELEMENT_H
#define FEATKHEX8ELEMENT_H

#include <featk/geometry/featkElement.h>

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension> class featkElement;

using featkHex8Element = featkElement<3, 8, 6, 3>;

#endif // FEATKHEX8ELEMENT_H
