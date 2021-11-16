/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featktIntegrationRuleBase.h

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
 * @class featkIntegrationRuleBase
 *
 * @brief Base class for numerical integration rules.
 *
 * featkIntegrationRuleBase is a base class for numerical integration rules.
 * featkIntegrationRuleBase is templated over unsigned integers NaturalCoordinates and
 * Points (see definitions below).
 *
 * Derived class constructor specializations must initialize the inheritted member
 * variable featkIntegrationRuleInterface::pointsAndWeights accordingly.
 *
 * @tparam NaturalDimension The natural dimension of the rule integration
 * points.
 *
 * @tparam Points The number of integration points of the rule.
 *
 */

#ifndef FEATKINTERGRATIONRULEBASE_H
#define FEATKINTERGRATIONRULEBASE_H

#include <featk/integration/featkIntegrationRuleInterface.h>

template<unsigned int NaturalDimension, unsigned int Points>
class featkIntegrationRuleBase : public featkIntegrationRuleInterface<NaturalDimension> {

    public:

        virtual ~featkIntegrationRuleBase();

    protected:

        featkIntegrationRuleBase();
};

template<unsigned int NaturalDimension, unsigned int Points>
featkIntegrationRuleBase<NaturalDimension, Points>::featkIntegrationRuleBase() {

}

template<unsigned int NaturalDimension, unsigned int Points>
featkIntegrationRuleBase<NaturalDimension, Points>::~featkIntegrationRuleBase() {

}

#endif // FEATKINTERGRATIONRULEBASE_H
