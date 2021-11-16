/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkNode.h

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
 * @class featkNode
 *
 * @brief Node with its attributes including coordinates in Dimension
 * cartesian dimensions.
 *
 * featkNode is the implementation of a node with its attributes including
 * coordinates in Dimension cartesian dimensions.
 *
 * featkNode stores a pointer to each featkElementInterface it belongs to
 * for computational efficiency when evaluating node derivative quantities
 * that need to be averaged over elements sharing this node.
 *
 * Each featkNode must be given a unique id at construction time whose
 * uniqueness is ensured by the user.
 *
 * @warning featkNode cartesian coordinates are automatically registered
 * as a node attribute in the featkMesh constructor under name "Cartesian
 * Coordinates". Make sure not to assign node attribute of order 1 with
 * name "Cartesian Coordinates" to featkMesh objects except if you intend
 * to modify the featkMesh geometry.
 *
 * @tparam Dimension The cartesian dimension of the node.
 *
 */

#ifndef FEATKNODE_H
#define FEATKNODE_H

#include <featk/core/featkDefines.h>
#include <featk/geometry/featkAttributable.h>
#include <featk/geometry/featkElementInterface.h>

#include <vector>

template<unsigned int Dimension> class featkElementInterface;

template<unsigned int Dimension>
class featkNode : public featkAttributable<Dimension> {

    public:

        featkNode(size_t id, AttributeValueType<Dimension, 1> coordinates);
        ~featkNode();

        void addElement(featkElementInterface<Dimension>* element);
        AttributeValueType<Dimension, 1> getCoordinates() const;
        std::vector<featkElementInterface<Dimension>*> getElements() const;
        size_t getID() const;

    private:

        std::vector<featkElementInterface<Dimension>*> elements;
        size_t id;  // Auto-assign id using a static counter variable?
};

template<unsigned int Dimension>
featkNode<Dimension>::featkNode(size_t id, AttributeValueType<Dimension, 1> coordinates) {

    this->id = id;
    this->setAttribute(1, std::make_shared<MatrixXd>(coordinates));
}

template<unsigned int Dimension>
featkNode<Dimension>::~featkNode() {

}


template<unsigned int Dimension>
void featkNode<Dimension>::addElement(featkElementInterface<Dimension>* element) {

    if (find(this->elements.begin(), this->elements.end(), element) == this->elements.end()) {

        this->elements.push_back(element);
    }
}

template<unsigned int Dimension>
AttributeValueType<Dimension, 1> featkNode<Dimension>::getCoordinates() const {

    return this->getAttributeValue(1);
}

template<unsigned int Dimension>
std::vector<featkElementInterface<Dimension>*> featkNode<Dimension>::getElements() const {

    return this->elements;
}

template<unsigned int Dimension>
size_t featkNode<Dimension>::getID() const {

    return this->id;
}

#endif // FEATKNODE_H
