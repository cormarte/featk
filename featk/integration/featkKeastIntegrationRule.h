/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featktKeastIntegrationRule.h

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
 * @class featkKeastIntegrationRule
 *
 * @brief Keast rule for numerical integration over the unit tetrahedron.
 *
 * featkKeastIntegrationRuleBase implements Keast rule for numerical
 * integration over a unit tetrahedron as proposed in P. Keast,
 * "Moderate-Degree Tetrahedral Quadrature Formulas," Comput. Methods. Appl.
 * Mech. Eng., vol. 55, no. 3, pp. 339-348, May 1986.
 * Constructor of featkKeastIntegrationRuleBase specializations initialize
 * the inheritted member variable
 * featkIntegrationRuleInterface::pointsAndWeights according to their Points
 * template parameter value.
 *
 * See J. E. Akin, Finite Element Analysis with Error Estimators: Chapter 10
 * - Integration Methods, p. 272, Table 10.4, 2005.
 *
 * @tparam Points The number of integration points of the rule.
 *
 */

#ifndef FEATKKEASTINTEGRATIONRULE_H
#define FEATKKEASTINTEGRATIONRULE_H

#include <featk/integration/featkIntegrationRuleBase.h>

template<unsigned int Points>
class featkKeastIntegrationRule final : public featkIntegrationRuleBase<3, Points> {

    public:

        featkKeastIntegrationRule();
        ~featkKeastIntegrationRule();
};

template<unsigned int Points>
featkKeastIntegrationRule<Points>::~featkKeastIntegrationRule() {

}

#endif // FEATKKEASTINTEGRATIONRULE_H
