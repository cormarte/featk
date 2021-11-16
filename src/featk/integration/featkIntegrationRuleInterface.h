/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featktIntegrationRuleInterface.h

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
 * @class featkIntegrationRuleInterface
 *
 * @brief Interface for featkIntegrationRuleBase objects with
 * NaturalDimensions natural dimensions.
 *
 * featkIntegrationRuleInterface is an interface allowing to handle
 * featkIntegrationRuleBase objects with NaturalDimension natural dimensions
 * and various Points template parameter values
 *
 * @tparam NaturalDimension The natural dimension of the integration points
 * of the rule.
 *
 * @tparam Points The number of integration points of the rule.
 *
 */

#ifndef FEATKINTEGRATIONRULEINTERFACE_H
#define FEATKINTEGRATIONRULEINTERFACE_H

#include <featk/core/featkDefines.h>

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

template<unsigned int NaturalDimension>
class featkIntegrationRuleInterface {

    public:

        virtual ~featkIntegrationRuleInterface();

        std::vector<std::pair<double, NaturalCoordinatesMatrixType<NaturalDimension>>> getPointsAndWeights();

    protected:

        featkIntegrationRuleInterface();

        std::vector<std::pair<double, NaturalCoordinatesMatrixType<NaturalDimension>>> pointsAndWeights;
};

template<unsigned int NaturalDimension>
featkIntegrationRuleInterface<NaturalDimension>::featkIntegrationRuleInterface() {

}

template<unsigned int NaturalDimension>
featkIntegrationRuleInterface<NaturalDimension>::~featkIntegrationRuleInterface() {

}

template<unsigned int NaturalDimension>
std::vector<std::pair<double, NaturalCoordinatesMatrixType<NaturalDimension>>> featkIntegrationRuleInterface<NaturalDimension>::getPointsAndWeights() {

    return this->pointsAndWeights;
}

#endif // FEATKINTEGRATIONRULEINTERFACE_H
