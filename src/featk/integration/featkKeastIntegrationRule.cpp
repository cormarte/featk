#include <featk/integration/featkKeastIntegrationRule.h>

featkKeastIntegrationRule<1>::featkKeastIntegrationRule() {

    this->pointsAndWeights = {std::make_pair(1.0/6.0, (NaturalCoordinatesMatrixType<3>() << 1.0/4.0, 1.0/4.0, 1.0/4.0).finished())};
}

featkKeastIntegrationRule<4>::featkKeastIntegrationRule() {

    this->pointsAndWeights = {std::make_pair(1.0/24.0, (NaturalCoordinatesMatrixType<3>() << (5.0+3.0*sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0).finished()),
                              std::make_pair(1.0/24.0, (NaturalCoordinatesMatrixType<3>() <<     (5.0-sqrt(5.0))/20.0, (5.0+3.0*sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0).finished()),
                              std::make_pair(1.0/24.0, (NaturalCoordinatesMatrixType<3>() <<     (5.0-sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0, (5.0+3.0*sqrt(5.0))/20.0).finished()),
                              std::make_pair(1.0/24.0, (NaturalCoordinatesMatrixType<3>() <<     (5.0-sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0,     (5.0-sqrt(5.0))/20.0).finished())};
}

featkKeastIntegrationRule<5>::featkKeastIntegrationRule() {

    this->pointsAndWeights = {std::make_pair(-4.0/30.0, (NaturalCoordinatesMatrixType<3>() << 1.0/4.0, 1.0/4.0, 1.0/4.0).finished()),
                              std::make_pair(9.0/120.0, (NaturalCoordinatesMatrixType<3>() << 1.0/2.0, 1.0/6.0, 1.0/6.0).finished()),
                              std::make_pair(9.0/120.0, (NaturalCoordinatesMatrixType<3>() << 1.0/6.0, 1.0/2.0, 1.0/6.0).finished()),
                              std::make_pair(9.0/120.0, (NaturalCoordinatesMatrixType<3>() << 1.0/6.0, 1.0/6.0, 1.0/2.0).finished()),
                              std::make_pair(9.0/120.0, (NaturalCoordinatesMatrixType<3>() << 1.0/6.0, 1.0/6.0, 1.0/6.0).finished())};
}

featkKeastIntegrationRule<11>::featkKeastIntegrationRule() {

    this->pointsAndWeights = {std::make_pair( -74.0/5625.0, (NaturalCoordinatesMatrixType<3>() <<                  1.0/4.0,                  1.0/4.0,                  1.0/4.0).finished()),
                              std::make_pair(343.0/45000.0, (NaturalCoordinatesMatrixType<3>() <<                11.0/14.0,                 1.0/14.0,                 1.0/14.0).finished()),
                              std::make_pair(343.0/45000.0, (NaturalCoordinatesMatrixType<3>() <<                 1.0/14.0,                11.0/14.0,                 1.0/14.0).finished()),
                              std::make_pair(343.0/45000.0, (NaturalCoordinatesMatrixType<3>() <<                 1.0/14.0,                 1.0/14.0,                11.0/14.0).finished()),
                              std::make_pair(343.0/45000.0, (NaturalCoordinatesMatrixType<3>() <<                 1.0/14.0,                 1.0/14.0,                 1.0/14.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0+sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0+sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0+sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0-sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0-sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0).finished()),
                              std::make_pair(  56.0/2250.0, (NaturalCoordinatesMatrixType<3>() << (1.0-sqrt(5.0/14.0))/4.0, (1.0-sqrt(5.0/14.0))/4.0, (1.0+sqrt(5.0/14.0))/4.0).finished())};
}
