/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkElementInterface.h

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
 * @class featkElementInterface
 *
 * @brief Interface for featkElement objects with Dimension cartesian
 * dimensions.
 *
 * featkElementInterface is an interface allowing to handle featkElement
 * objects with Dimension cartesian dimensions and various Nodes, Boundaries
 * and NaturalDimension template parameter values. featkElementInterface
 * makes it possible to define hybrid featkMesh seamlessly.
 *
 * @tparam Dimension The cartesian dimension of the element.
 *
 */

#ifndef FEATKELEMENTINTERFACE_H
#define FEATKELEMENTINTERFACE_H

#include <featk/core/featkDefines.h>
#include <featk/geometry/featkAttributable.h>
#include <featk/geometry/featkHex8Element.h>
#include <featk/geometry/featkNode.h>
//#include <featk/geometry/featkTet4Element.h>

#include <vector>

template<unsigned int Dimension> class featkNode;

template<unsigned int Dimension>
class featkElementInterface : public featkAttributable<Dimension> {

    public:

        virtual ~featkElementInterface();

        template<unsigned int Order> MatrixXd getBtBIntegralMatrix() const;
        template<unsigned int Order> MatrixXd getBtCBIntegralMatrix(size_t elementAttributeID) const;
        template<unsigned int Order> VectorXd getNodeBQVector(featkNode<Dimension>* node, size_t elementAttributeID) const;
        template<unsigned int Order> MatrixXd getNtCNIntegralMatrix(size_t elementAttributeID) const;
        template<unsigned int Order> VectorXd getNtCNQNQIntegralVector(size_t elementAttributeID, size_t nodeAttributeID) const;
        template<unsigned int Order> MatrixXd getNtNIntegralMatrix() const;
        template<unsigned int Order> VectorXd getNtNQIntegralVector(size_t nodeAttributeID) const;
        template<unsigned int Order> VectorXd getNtNQNQIntegralVector(size_t nodeAttributeID) const;
        template<unsigned int Order> VectorXd getQVector(size_t nodeAttributeID) const;

        AttributeValueType<Dimension, 1> getBarycenter() const;
        featkElementType getElementType() const;
        featkNode<Dimension>* getNode(unsigned int index) const;
        std::vector<featkNode<Dimension>*> getNodes() const;

    protected:

        featkElementInterface(featkElementType type);

        const featkElementType elementType;
        std::vector<featkNode<Dimension>*> nodes;                                               // Or use std::Array<featkNode<Dimension>, Nodes> in featkElementBase and virtual std::vector getNodes featkElementInterface?
};

template<unsigned int Dimension>
featkElementInterface<Dimension>::featkElementInterface(featkElementType type) : elementType(type) {

}

template<unsigned int Dimension>
featkElementInterface<Dimension>::~featkElementInterface() {

}


template<unsigned int Dimension>
template<unsigned int Order>
MatrixXd featkElementInterface<Dimension>::getBtBIntegralMatrix() const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getBtBIntegralMatrix<Order>();
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getBtBIntegralMatrix<Order>();
            break;

        default:
            matrix = MatrixXd::Zero(1, 1);;
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
MatrixXd featkElementInterface<Dimension>::getBtCBIntegralMatrix(size_t elementAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getBtCBIntegralMatrix<Order>(elementAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getBtCBIntegralMatrix<Order>(elementAttributeID);
            break;

        default:
            matrix = MatrixXd::Zero(1, 1);;
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
VectorXd featkElementInterface<Dimension>::getNodeBQVector(featkNode<Dimension>* node, size_t elementAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getNodeBQMatrix<Order>(node, elementAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getNodeBQMatrix<Order>(node, elementAttributeID);
            break;

        default:
            matrix = VectorXd::Zero(1, 1);;
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
MatrixXd featkElementInterface<Dimension>::getNtCNIntegralMatrix(size_t elementAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getNtCNIntegralMatrix<Order>(elementAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getNtCNIntegralMatrix<Order>(elementAttributeID);
            break;

        default:
            matrix = MatrixXd::Zero(1, 1);;
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
MatrixXd featkElementInterface<Dimension>::getNtNIntegralMatrix() const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getNtNIntegralMatrix<Order>();
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getNtNIntegralMatrix<Order>();
            break;

        default:
            matrix = MatrixXd::Zero(1, 1);
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
VectorXd featkElementInterface<Dimension>::getNtNQIntegralVector(size_t nodeAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getNtNQIntegralMatrix<Order>(nodeAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getNtNQIntegralMatrix<Order>(nodeAttributeID);
            break;


        default:
            matrix = VectorXd::Zero(1);
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
VectorXd featkElementInterface<Dimension>::getNtCNQNQIntegralVector(size_t elementAttributeID, size_t nodeAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getNtCNQNQIntegralMatrix<Order>(elementAttributeID, nodeAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getNtCNQNQIntegralMatrix<Order>(elementAttributeID, nodeAttributeID);
            break;


        default:
            matrix = VectorXd::Zero(1);
            break;
    }

    return matrix;
}

template<unsigned int Dimension>
template<unsigned int Order>
VectorXd featkElementInterface<Dimension>::getQVector(size_t nodeAttributeID) const {

    /**
     * C++ does not allow template virtual member function but as many matrix dimensions as possible must be
     * known at compile time for performance reasons, hence element matrix getters remain templated over
     * variable order. The following solution is not elegant but does the trick while avoiding multiple
     * dynamic_cast attempts.
     *
     * Another solution would be to not allow hybrid meshes and template featkMesh over element type and
     * featkSolverBase over featkMesh<Dimension, Element>. This way, all matrix dimension would be known at
     * compile time.
     */

    MatrixXd matrix;

    switch (this->elementType) {

        case FEATK_TET4:
            matrix = static_cast<const featkTet4Element*>(this)->getQMatrix<Order>(nodeAttributeID);
            break;

        case FEATK_HEX8:
            matrix = static_cast<const featkHex8Element*>(this)->getQMatrix<Order>(nodeAttributeID);
            break;

        default:
            matrix = VectorXd::Zero(1);
            break;
    }

    return matrix;
}


template<unsigned int Dimension>
AttributeValueType<Dimension, 1> featkElementInterface<Dimension>::getBarycenter() const {

    AttributeValueType<Dimension, 1> barycenter = AttributeValueType<Dimension, 1>::Zero();

    for (featkNode<Dimension>* node : this->nodes) {

        barycenter += node->getCoordinates();
    }

    barycenter /= n;

    return barycenter;
}

template<unsigned int Dimension>
featkElementType featkElementInterface<Dimension>::getElementType() const {

    return this->elementType;
}

template<unsigned int Dimension>
featkNode<Dimension>* featkElementInterface<Dimension>::getNode(unsigned int index) const {

    //return (index < this->getNumberOfNodes()) ? this->nodes[index] : nullptr;
    return (index < this->nodes.size()) ? this->nodes[index] : nullptr;
}

template<unsigned int Dimension>
std::vector<featkNode<Dimension>*> featkElementInterface<Dimension>::getNodes() const {

    return this->nodes;
}

#endif // FEATKELEMENTINTERFACE_H
