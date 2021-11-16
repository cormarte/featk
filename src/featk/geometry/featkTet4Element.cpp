#include <featk/geometry/featkTet4Element.h>
#include <featk/integration/featkKeastIntegrationRule.h>

const featkTet4Element::NodesNaturalCoordinatesMatrixType featkTet4Element::nodesNaturalCoordinates = (featkTet4Element::NodesNaturalCoordinatesMatrixType() << 0.0, 0.0, 0.0,
                                                                                                                                                                1.0, 0.0, 0.0,
                                                                                                                                                                0.0, 1.0, 0.0,
                                                                                                                                                                0.0, 0.0, 1.0).finished();

featkIntegrationRuleInterface<3>* featkTet4Element::integrationRule = new featkKeastIntegrationRule<4>();

featkTet4Element::featkElement(std::vector<featkNode<3>*> nodes) : featkElementInterface<3>(FEATK_TET4) {

    this->nodes = nodes;  // Check number of nodes

    for (unsigned int i=0; i!=4; i++) {

        this->nodes[i]->addElement(this);
    }
}

featkTet4Element::ShapeFunctionNaturalDerivativeValuesMatrixType featkTet4Element::getShapeFunctionNaturalDerivativeValues(const NaturalCoordinatesMatrixType& point) const {

    return (featkTet4Element::ShapeFunctionNaturalDerivativeValuesMatrixType() << -1.0, 1.0, 0.0, 0.0,
                                                                                  -1.0, 0.0, 1.0, 0.0,
                                                                                  -1.0, 0.0, 0.0, 1.0).finished();
}

featkTet4Element::ShapeFunctionValuesMatrixType featkTet4Element::getShapeFunctionValues(const NaturalCoordinatesMatrixType& point) const {

    return (ShapeFunctionValuesMatrixType() << 1.0-point(0, 0)-point(0, 1)-point(0, 2), point(0, 0), point(0,1), point(0, 2)).finished();
}
