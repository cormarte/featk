#include <featk/geometry/featkHex8Element.h>
#include <featk/integration/featkProductGaussianQuadratureIntegrationRule.h>

const featkHex8Element::NodesNaturalCoordinatesMatrixType featkHex8Element::nodesNaturalCoordinates = (featkHex8Element::NodesNaturalCoordinatesMatrixType() <<  -1.0, -1.0, -1.0,
                                                                                                                                                                  1.0, -1.0, -1.0,
                                                                                                                                                                  1.0,  1.0, -1.0,
                                                                                                                                                                 -1.0,  1.0, -1.0,
                                                                                                                                                                 -1.0, -1.0,  1.0,
                                                                                                                                                                  1.0, -1.0,  1.0,
                                                                                                                                                                  1.0,  1.0,  1.0,
                                                                                                                                                                 -1.0,  1.0,  1.0).finished();

featkIntegrationRuleInterface<3>* featkHex8Element::integrationRule = new featkProductGaussianQuadratureIntegrationRule<3, 8>();

featkHex8Element::featkElement(std::vector<featkNode<3>*> nodes) : featkElementInterface<3>(FEATK_HEX8) {

    this->nodes = nodes;  // Check number of nodes

    for (unsigned int i=0; i!=8; i++) {

        this->nodes[i]->addElement(this);
    }
}

featkHex8Element::ShapeFunctionNaturalDerivativeValuesMatrixType featkHex8Element::getShapeFunctionNaturalDerivativeValues(const NaturalCoordinatesMatrixType& point) const {

    ShapeFunctionNaturalDerivativeValuesMatrixType derivatives;

    for (unsigned int n=0; n!=8; n++) {

        derivatives(0, n) = (1.0/8.0)*featkHex8Element::nodesNaturalCoordinates(n,0)*(1.0+featkHex8Element::nodesNaturalCoordinates(n,1)*point(0,1))*(1.0+featkHex8Element::nodesNaturalCoordinates(n,2)*point(0,2));
        derivatives(1, n) = (1.0/8.0)*featkHex8Element::nodesNaturalCoordinates(n,1)*(1.0+featkHex8Element::nodesNaturalCoordinates(n,0)*point(0,0))*(1.0+featkHex8Element::nodesNaturalCoordinates(n,2)*point(0,2));
        derivatives(2, n) = (1.0/8.0)*featkHex8Element::nodesNaturalCoordinates(n,2)*(1.0+featkHex8Element::nodesNaturalCoordinates(n,0)*point(0,0))*(1.0+featkHex8Element::nodesNaturalCoordinates(n,1)*point(0,1));
    }

    return derivatives;
}

featkHex8Element::ShapeFunctionValuesMatrixType featkHex8Element::getShapeFunctionValues(const NaturalCoordinatesMatrixType& point) const {

    ShapeFunctionValuesMatrixType values;

    for (unsigned int n=0; n!=8; n++) {

        values(0, n) = (1.0/8.0)*(1.0+featkHex8Element::nodesNaturalCoordinates(n, 0)*point(0, 0))*(1.0+featkHex8Element::nodesNaturalCoordinates(n, 1)*point(0, 1))*(1.0+featkHex8Element::nodesNaturalCoordinates(n, 2)*point(0, 2));
    }

    return values;
}
