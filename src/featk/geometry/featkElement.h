/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkElement.h

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
 * @class featkElement
 *
 * @brief Concrete implementation of a featkElementInterface.
 *
 * featkElement is a the concrete implementation of an element with its
 * attributes, nodes, shape function values and derivatives and
 * integration rule. featkElement provides routine member functions for
 * numerical integration of usual quantities found in finite element
 * problems.
 *
 * featkElement is templated over unsigned integers Dimension, Nodes,
 * Boundaries and NaturalDimension (see definitions below). A given
 * combination of these template parameters univoquely defines an
 * element type. Member variables featkElement::integrationRule and
 * featkElement::nodesNaturalCoordinates and member functions
 * featkElement::featkElement(),
 * featkElement::getShapeFunctionNaturalDerivativeValues() and
 * featkElement::getShapeFunctionValues() have no default implementation and
 * must be reimplemented for each template specialization.
 *
 * @tparam Dimension The cartesian dimension of the element.
 *
 * @tparam Nodes The number of nodes of the element.
 *
 * @tparam Boundaries The number of boundaries of the element.
 *
 * @tparam NaturalDimension The natural dimension of the element
 * parametrization.
 *
 */

#ifndef FEATKELEMENT_H
#define FEATKELEMENT_H

#include <featk/core/featkDefines.h>
#include <featk/core/featkGlobal.h>
#include <featk/geometry/featkElementInterface.h>
#include <featk/integration/featkIntegrationRuleInterface.h>

#include <vector>

template<unsigned int Dimension> class featkElementInterface;
template<unsigned int Dimension> class featkNode;

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
class featkElement final : public featkElementInterface<Dimension> {

    public:

        template<unsigned int Order> using BMatrixType = BMatrixType<Dimension, Nodes, Order>;
        template<unsigned int Order> using CMatrixType = CMatrixType<Dimension, Order>;
        template<unsigned int Order> using DMatrixType = DMatrixType<Dimension, Order>;
        template<unsigned int Order> using KMatrixType = KMatrixType<Dimension, Nodes, Order>;
        template<unsigned int Order> using NMatrixType = NMatrixType<Dimension, Nodes, Order>;
        template<unsigned int Order> using QMatrixType = QMatrixType<Dimension, Nodes, Order>;

        using JacobianMatrixType                               = JacobianMatrixType<Dimension, NaturalDimension>;
        using NaturalCoordinatesMatrixType                     = NaturalCoordinatesMatrixType<NaturalDimension>;
        using NodesCartesianCoordinatesMatrixType              = NodesCartesianCoordinatesMatrixType<Dimension, Nodes>;
        using NodesNaturalCoordinatesMatrixType                = NodesNaturalCoordinatesMatrixType <Nodes, NaturalDimension>;
        using ShapeFunctionCartesianDerivativeValuesMatrixType = ShapeFunctionCartesianDerivativeValuesMatrixType<Dimension, Nodes>;
        using ShapeFunctionNaturalDerivativeValuesMatrixType   = ShapeFunctionNaturalDerivativeValuesMatrixType<Nodes, NaturalDimension>;
        using ShapeFunctionValuesMatrixType                    = ShapeFunctionValuesMatrixType<Nodes>;

        FEATK_EXPORT featkElement(std::vector<featkNode<Dimension>*> nodes);
        ~featkElement();

        FEATK_EXPORT ShapeFunctionNaturalDerivativeValuesMatrixType getShapeFunctionNaturalDerivativeValues(const NaturalCoordinatesMatrixType& point) const;  // May be known at compile time for each integration point if integration rule is given as template parameters
        FEATK_EXPORT ShapeFunctionValuesMatrixType getShapeFunctionValues(const NaturalCoordinatesMatrixType& point) const;                                    // May be known at compile time for each integration point if integration rule is given as template parameters

        static void setIntegrationRule(featkIntegrationRuleInterface<Dimension>* rule);
        static featkIntegrationRuleInterface<Dimension> getIntegrationRule();
        static NodesNaturalCoordinatesMatrixType getNodesNaturalCoordinates();

        template<unsigned int Order> BMatrixType<Order> getBMatrix(const ShapeFunctionCartesianDerivativeValuesMatrixType& cartesianDerivatives) const;
        template<unsigned int Order> KMatrixType<Order> getBtBIntegralMatrix() const;
        template<unsigned int Order> KMatrixType<Order> getBtCBIntegralMatrix(size_t elementAttributeID) const;
        template<unsigned int Order> CMatrixType<Order> getCMatrix(size_t elementAttributeID) const;
        template<unsigned int Order> NMatrixType<Order> getNMatrix(const NaturalCoordinatesMatrixType& point) const;
        template<unsigned int Order> DMatrixType<Order> getNodeBQMatrix(featkNode<Dimension>* node, size_t nodeAttributeID) const;
        template<unsigned int Order> DMatrixType<Order> getNodeCBQMatrix(featkNode<Dimension>* node, size_t elementAttributeID, size_t nodeAttributeID) const;
        template<unsigned int Order> KMatrixType<Order> getNtCNIntegralMatrix(size_t elementAttributeID) const;
        template<unsigned int Order> KMatrixType<Order> getNtNIntegralMatrix() const;
        template<unsigned int Order> QMatrixType<Order> getNtNQIntegralMatrix(size_t nodeAttributeID) const;
        template<unsigned int Order> QMatrixType<Order> getNtCNQNQIntegralMatrix(size_t elementAttributeID, size_t nodeAttributeID) const;
        template<unsigned int Order> QMatrixType<Order> getQMatrix(size_t nodeAttributeID) const;

        JacobianMatrixType getJacobian(const NaturalCoordinatesMatrixType& point) const;
        JacobianMatrixType getJacobian(const ShapeFunctionNaturalDerivativeValuesMatrixType& naturalDerivatives) const;
        NodesCartesianCoordinatesMatrixType getNodeCartesianCoordinates() const;
        ShapeFunctionCartesianDerivativeValuesMatrixType getShapeFunctionCartesianDerivativeValues(const NaturalCoordinatesMatrixType& point) const;
        ShapeFunctionCartesianDerivativeValuesMatrixType getShapeFunctionCartesianDerivativeValues(const ShapeFunctionNaturalDerivativeValuesMatrixType& naturalDerivatives, const JacobianMatrixType& jacobian) const;

    private:

        FEATK_EXPORT static featkIntegrationRuleInterface<Dimension>* integrationRule;
        FEATK_EXPORT static const NodesNaturalCoordinatesMatrixType nodesNaturalCoordinates;
};

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::~featkElement() {

}


template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
void featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::setIntegrationRule(featkIntegrationRuleInterface<Dimension>* rule) {

    featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::integrationRule = rule;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
featkIntegrationRuleInterface<Dimension> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getIntegrationRule() {

    return featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::integrationRule;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::NodesNaturalCoordinatesMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNodesNaturalCoordinates() {

    return featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::nodesNaturalCoordinates;
}


template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::BMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getBMatrix(const ShapeFunctionCartesianDerivativeValuesMatrixType& cartesianDerivatives) const {

    /**
     *                                                        T                                        T
     * [u_x,x u_x,y u_x,z u_y,x u_y,y u_y,z u_z,x u_z,y u_z,z]  = B*[u_1x u_1y u_1z ... u_nx u_ny u_nz]
     */

    BMatrixType<Order> b = BMatrixType<Order>::Zero();

    for (int i=0; i!=POWER(Dimension, Order+1); i++) {

        for (int j=0; j!=Nodes; j++) {

            b(i, j*POWER(Dimension, Order)+i/Dimension) = cartesianDerivatives(i%Dimension, j);
        }
    }

    return b;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::KMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getBtBIntegralMatrix() const {

    KMatrixType<Order> k = KMatrixType<Order>::Zero();

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(point);
        JacobianMatrixType jacobian = this->getJacobian(naturalDerivatives);
        ShapeFunctionCartesianDerivativeValuesMatrixType cartesianDerivatives = this->getShapeFunctionCartesianDerivativeValues(naturalDerivatives, jacobian);
        BMatrixType<Order> b = this->getBMatrix<Order>(cartesianDerivatives);

        k += weight*jacobian.determinant()*b.transpose()*b;
    }

    return k;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::KMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getBtCBIntegralMatrix(size_t elementAttributeID) const {

    KMatrixType<Order> k = KMatrixType<Order>::Zero();
    CMatrixType<2*(Order+1)> c = this->getCMatrix<2*(Order+1)>(elementAttributeID);

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        /**
          * In the previous implementation, Jacobian (and shape function natural derivatives) were computed twice.
          */

        /*
        NaturalCoordinatesMatrixType point = integrationPoints[p];
        BMatrixType<Order> b = this->getBMatrix<Order>(point);
        double jacobianDeterminant = this->getJacobian(point).determinant();

         k += weights[p]*jacobianDeterminant*b.transpose()*(c*b);
        */

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(point);
        JacobianMatrixType jacobian = this->getJacobian(naturalDerivatives);
        ShapeFunctionCartesianDerivativeValuesMatrixType cartesianDerivatives = this->getShapeFunctionCartesianDerivativeValues(naturalDerivatives, jacobian);
        BMatrixType<Order> b = this->getBMatrix<Order>(cartesianDerivatives);

        k += weight*jacobian.determinant()*b.transpose()*c*b;
    }

    return k;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::CMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getCMatrix(size_t elementAttributeID) const {

    CMatrixType<Order> c = this->getAttributeValue(elementAttributeID);

    return c;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::NMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNMatrix(const NaturalCoordinatesMatrixType& point) const {

    /**
     *              T                                        T
     * [u_x u_y u_z]  = N*[u_1x u_1y u_1z ... u_nx u_ny u_nz]
     */

    NMatrixType<Order> n = NMatrixType<Order>::Zero();
    ShapeFunctionValuesMatrixType values = this->getShapeFunctionValues(point);

    for (int i=0; i!=POWER(Dimension, Order); i++) {

        for (int j=0; j!=Nodes; j++) {

            n(i, i+j*POWER(Dimension, Order)) = values(0, j);
        }
    }

    return n;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::DMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNodeBQMatrix(featkNode<Dimension>* node, size_t nodeAttributeID) const {

    size_t nodeIndex = distance(this->nodes.begin(), find(this->nodes.begin(), this->nodes.end(), node));

    ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(this->nodesNaturalCoordinates.row(nodeIndex));
    JacobianMatrixType jacobian = this->getJacobian(naturalDerivatives);
    ShapeFunctionCartesianDerivativeValuesMatrixType cartesianDerivatives = this->getShapeFunctionCartesianDerivativeValues(naturalDerivatives, jacobian);
    BMatrixType<Order> b = this->getBMatrix<Order>(cartesianDerivatives);

    QMatrixType<Order> q = this->getQMatrix<Order>(nodeAttributeID);

    return b*q;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::DMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNodeCBQMatrix(featkNode<Dimension>* node, size_t elementAttributeID, size_t nodeAttributeID) const {

    size_t nodeIndex = distance(this->nodes.begin(), find(this->nodes.begin(), this->nodes.end(), node));

    CMatrixType<2*(Order+1)> c = this->getCMatrix<2*(Order+1)>(elementAttributeID);

    ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(this->nodesNaturalCoordinates.row(nodeIndex));
    JacobianMatrixType jacobian = this->getJacobian(naturalDerivatives);
    ShapeFunctionCartesianDerivativeValuesMatrixType cartesianDerivatives = this->getShapeFunctionCartesianDerivativeValues(naturalDerivatives, jacobian);
    BMatrixType<Order> b = this->getBMatrix<Order>(cartesianDerivatives);

    QMatrixType<Order> q = this->getQMatrix<Order>(nodeAttributeID);

    return c*b*q;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::KMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNtCNIntegralMatrix(size_t elementAttributeID) const {

    KMatrixType<Order> k = KMatrixType<Order>::Zero();
    CMatrixType<Order> c = this->getCMatrix<Order>(elementAttributeID);

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        NMatrixType<Order> n = this->getNMatrix<Order>(point);
        JacobianMatrixType jacobian = this->getJacobian(point);

        k += weight*jacobian.determinant()*n.transpose()*c*n;
    }

    return k;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::KMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNtNIntegralMatrix() const {

    KMatrixType<Order> k = KMatrixType<Order>::Zero();

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        NMatrixType<Order> n = this->getNMatrix<Order>(point);
        JacobianMatrixType jacobian = this->getJacobian(point);

        k += weight*jacobian.determinant()*n.transpose()*n;
    }

    return k;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::QMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNtNQIntegralMatrix(size_t nodeAttributeID) const {

    QMatrixType<Order> f = QMatrixType<Order>::Zero();
    QMatrixType<Order> q = this->getQMatrix<Order>(nodeAttributeID);

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        NMatrixType<Order> n = this->getNMatrix<Order>(point);
        JacobianMatrixType jacobian = this->getJacobian(point);

        f += weight*jacobian.determinant()*n.transpose()*n*q;
    }

    return f;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::QMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNtCNQNQIntegralMatrix(size_t elementAttributeID, size_t nodeAttributeID) const {

    QMatrixType<Order> f = QMatrixType<Order>::Zero();
    CMatrixType<Order> c = this->getCMatrix<Order>(elementAttributeID);
    QMatrixType<Order> q = this->getQMatrix<Order>(nodeAttributeID);

    std::vector<std::pair<double, NaturalCoordinatesMatrixType>> pointsAndWeights = this->integrationRule->getPointsAndWeights();

    for (std::pair<double, NaturalCoordinatesMatrixType> p : pointsAndWeights) {

        double weight = p.first;
        NaturalCoordinatesMatrixType point = p.second;

        NMatrixType<Order> n = this->getNMatrix<Order>(point);
        JacobianMatrixType jacobian = this->getJacobian(point);

        f += weight*jacobian.determinant()*n.transpose()*c*n*q*n*q;
    }

    return f;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
template<unsigned int Order>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::QMatrixType<Order> featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getQMatrix(size_t nodeAttributeID) const {

    /**
     * The following implementation should be more efficient, however Eigen does not allow to define Matrix<double, row, 1, RowMajor>
     * types (i.e featkAttribute<Dimension, 1>::ValueType types with row major storage), hence ColMajor matrices must be used for
     * featkAttribute<Dimension, Order>::ValueType types which are reshaped to Matrix<double, POWER(Dimension, Order), 1> column
     * by column by Eigen::Map, i.e.:
     *
     *           |a|
     *  |a b|    |c|
     *  |c d| => |b|
     *           |d|
     *
     * This is however not a problem for symmetric attributes.
     */

    /*
    QMatrixType<Order> q;

    for (size_t n=0; n!=Nodes; n++) {

        q.block(POWER(Dimension, Order)*n, 0, POWER(Dimension, Order), 1) = Map<Matrix<double, POWER(Dimension, Order), 1>>(this->nodes[n]->getAttributeValue<Order>(nodeAttributeName).data());
    }

    return q;
    */

    QMatrixType<Order> q;
    size_t a = 0;

    for (unsigned int n=0; n!=Nodes; n++) {

        AttributeValueType<Dimension, Order> value = this->nodes[n]->getAttributeValue(nodeAttributeID);

        for (size_t i=0; i!=AttributeValueType<Dimension, Order>::RowsAtCompileTime; i++) {

            q.block(a, 0, AttributeValueType<Dimension, Order>::ColsAtCompileTime, 1) = value.row(i).transpose();
            a += AttributeValueType<Dimension, Order>::ColsAtCompileTime;
        }
    }

    return q;
}


template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::JacobianMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getJacobian(const NaturalCoordinatesMatrixType& point) const {

    ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(point);

    return this->getJacobian(naturalDerivatives);
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::JacobianMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getJacobian(const ShapeFunctionNaturalDerivativeValuesMatrixType& naturalDerivatives) const {

    NodesCartesianCoordinatesMatrixType cartesianCoordinates = this->getNodeCartesianCoordinates();
    JacobianMatrixType jacobian = (naturalDerivatives*cartesianCoordinates);

    return jacobian;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::NodesCartesianCoordinatesMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getNodeCartesianCoordinates() const {

    NodesCartesianCoordinatesMatrixType c;

    for (unsigned int n=0; n!=Nodes; n++) {

        c.row(n) = this->nodes[n]->getCoordinates().transpose();
    }

    return c;
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::ShapeFunctionCartesianDerivativeValuesMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getShapeFunctionCartesianDerivativeValues(const NaturalCoordinatesMatrixType& point) const {

    /**
     * See Felippa. 2017. Advanced Finite Element Methods. Ch.11 p.7. and Filomeno Coelho and Pyl. 2014. Structural Analysis and Finite Elements, p.71.
     */

    ShapeFunctionNaturalDerivativeValuesMatrixType naturalDerivatives = this->getShapeFunctionNaturalDerivativeValues(point);
    JacobianMatrixType jacobian = this->getJacobian(point);

    return this->getShapeFunctionCartesianDerivativeValues(naturalDerivatives, jacobian);
}

template<unsigned int Dimension, unsigned int Nodes, unsigned int Boundaries, unsigned int NaturalDimension=Dimension>
typename featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::ShapeFunctionCartesianDerivativeValuesMatrixType featkElement<Dimension, Nodes, Boundaries, NaturalDimension>::getShapeFunctionCartesianDerivativeValues(const ShapeFunctionNaturalDerivativeValuesMatrixType& naturalDerivatives, const JacobianMatrixType& jacobian) const {

    /**
     * See Felippa. 2017. Advanced Finite Element Methods. Ch.11 p.7. and Filomeno Coelho and Pyl. 2014. Structural Analysis and Finite Elements, p.71.
     */

    ShapeFunctionCartesianDerivativeValuesMatrixType cartesianDerivatives = jacobian.inverse()*naturalDerivatives;

    return cartesianDerivatives;
}

#endif // FEATKELEMENT_H
