#ifndef FEATKPRODUCTGAUSSIANQUADRATUREINTEGRATIONRULE_H
#define FEATKPRODUCTGAUSSIANQUADRATUREINTEGRATIONRULE_H

#include <featk/integration/featkIntegrationRuleBase.h>

#include <map>

template<unsigned int NaturalDimension, unsigned int Points>
class featkProductGaussianQuadratureIntegrationRule final : public featkIntegrationRuleBase<NaturalDimension, Points> {

    public:

        featkProductGaussianQuadratureIntegrationRule();
        ~featkProductGaussianQuadratureIntegrationRule();

    private:

        //static const std::map<unsigned int, std::vector<std::pair<double, double>>> pointsAndWeights1D;  // Prone to static initialization order fiasco
};

/*
template<unsigned int NaturalDimension, unsigned int Points>
const std::map<unsigned int, std::vector<std::pair<double, double>>> featkProductGaussianQuadratureIntegrationRule<NaturalDimension, Points>::pointsAndWeights1D = {{1, {{2.0, 0.0}}},

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
*/

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
