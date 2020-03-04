/*=========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkNode.h

  Copyright (c) Corentin Martens
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class featkNode
 * @brief element node with its attributes including cartesian coordinates.
 *
 * featkNode encapsulates element node attributes including its cartesian
 * coordinates in a Dimension dimensional space.
 *
 * Attributes are stored as Eigen::MatrixXd objects for ease but are
 * accessed through Eigen::Matrix objects whose dimensions are fixed at
 * compile time given the attribute order (0 for scalars, 1 for vectors, n
 * for nth-order tensor). This ensures the adequacy of most Eigen::Matrix
 * dimensions at compile time, limiting the risk of run time errors while
 * benefiting from the template construction of Eigen.
 *
 * featkNode object also store a pointer to each featkIElement it belongs
 * to for computational efficiency.
 *
 * Each node must be given a unique id at construction time whose uniqueness
 * is ensured by the user.
 *
 * @tparam Dimension the number of cartesian dimensions
 * of the node (1, 2, or 3).
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
