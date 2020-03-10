/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkProductGaussianQuadratureIntegrationRule.h

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
 * @class featkProductGaussianQuadratureIntegrationRule
 *
 * @brief Product Gaussian quadrature rule for numerical integration over a
 * reference element.
 *
 * featkProductGaussianQuadratureIntegrationRule implements product Gaussian
 * quadrature rule for numerical integration over a reference element in
 * natural coordinates. Integration points and weights are defined in
 * NaturalDimension 1 and are computed using the product rule for
 * NaturalDimension > 1.
 *
 * Constructor of featkProductGaussianQuadratureIntegrationRule
 * specializations initialize the inheritted member variable
 * featkIntegrationRuleInterface::pointsAndWeights according to their
 * template parameter values.
 *
 * See C. Felippa, Advanced Finite Element Methods: Chapter 11 - The 8-Node
 * Hexahedron, p. 10-11, 2017.
 *
 * @warning For a given NaturalDimension template parameter value, the
 * Points template parameter value can only be 1, 2, 3, 4 or 5 to the power
 * of NaturalDimension.
 *
 * @tparam NaturalDimension The natural dimension of the rule integration
 * points.
 *
 * @tparam Points The number of integration points of the rule.
 *
 */

#ifndef FEATKPRODUCTGAUSSIANQUADRATUREINTEGRATIONRULE_H
#define FEATKPRODUCTGAUSSIANQUADRATUREINTEGRATIONRULE_H

#include <featk/integration/featkIntegrationRuleBase.h>

#include <map>

template<unsigned int NaturalDimension, unsigned int Points>
class featkProductGaussianQuadratureIntegrationRule final : public featkIntegrationRuleBase<NaturalDimension, Points> {

    public:

        featkProductGaussianQuadratureIntegrationRule();
        ~featkProductGaussianQuadratureIntegrationRule();
};

template<unsigned int NaturalDimension, unsigned int Points>
featkProductGaussianQuadratureIntegrationRule<NaturalDimension, Points>::featkProductGaussianQuadratureIntegrationRule() {

    const std::map<unsigned int, std::vector<std::pair<double, double>>> pointsAndWeights1D = {{1, {{2.0, 0.0}}},

                                                                                               {2, {{1.0, -1.0/sqrt(3.0)},
                                                                                                    {1.0,  1.0/sqrt(3.0)}}},

                                                                                               {3, {{5.0/9.0, -sqrt(3.0/5.0)},
                                                                                                    {8.0/9.0,            0.0},
                                                                                                    {5.0/9.0,  sqrt(3.0/5.0)}}},

                                                                                               {4, {{(1.0/2.0)-sqrt(5.0/6.0)/6.0, -sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)},
                                                                                                    {(1.0/2.0)+sqrt(5.0/6.0)/6.0, -sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)},
                                                                                                    {(1.0/2.0)+sqrt(5.0/6.0)/6.0,  sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)},
                                                                                                    {(1.0/2.0)-sqrt(5.0/6.0)/6.0,  sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)}}},

                                                                                               {5, {{(322.0-13.0*sqrt(70.0))/900.0, -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0},
                                                                                                    {(322.0+13.0*sqrt(70.0))/900.0, -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0},
                                                                                                    {                  512.0/900.0,                               0.0},
                                                                                                    {(322.0+13.0*sqrt(70.0))/900.0,  sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0},
                                                                                                    {(322.0-13.0*sqrt(70.0))/900.0,  sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0}}}};

    unsigned int pointsPerDimension = (unsigned int)(std::pow(Points, 1.0/NaturalDimension));

    for (unsigned int p=0; p!=Points; p++) {

        NaturalCoordinatesMatrixType<NaturalDimension> point = NaturalCoordinatesMatrixType<NaturalDimension>();
        double weight = 1.0;

        unsigned int index = p;

        for (unsigned int d=0; d!=NaturalDimension; d++) {

            unsigned int n = POWER(pointsPerDimension, NaturalDimension-d-1);
            unsigned int i = index/n;

            weight *= pointsAndWeights1D.at(pointsPerDimension)[i].first;
            point(0, d) = pointsAndWeights1D.at(pointsPerDimension)[i].second;

            index -= i*n;
        }

        this->pointsAndWeights.push_back({weight, point});
    }
}

template<unsigned int NaturalDimension, unsigned int Points>
featkProductGaussianQuadratureIntegrationRule<NaturalDimension, Points>::~featkProductGaussianQuadratureIntegrationRule() {

}

#endif // FEATKPRODUCTGAUSSIANQUADRATUREINTEGRATIONRULE_H
