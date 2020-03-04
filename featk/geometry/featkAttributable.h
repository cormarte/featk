/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkAttributable.h

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
 * @class featkAttributable
 *
 * @brief Base class for featk objects possessing attributes, i.e. featkNode
 * and featkElementInterface.
 *
 * featkAttributable is a base class for featk objects with attributes
 * (featkNode and featkElementInterface), providing member variables and
 * functions for attribute assignment, storage and access.
 *
 * Attributes are stored as std::shared_ptr<Eigen::MatrixXd>.
 * std::shared_ptr allows featkAttributable objects to share attributes
 * without having to duplicate them nor to worry about memory mangement.
 * Dynamic-sized Eigen::MatrixXd allows to store attributes of various
 * Order seamlessly.
 *
 * @tparam Dimension The cartesian dimension of the attributable.
 *
 */

#ifndef FEATKATTRIBUTABLE_H
#define FEATKATTRIBUTABLE_H

#include <Eigen/Dense>
#include <map>
#include <memory>

using namespace Eigen;

template<unsigned int Dimension> class featkMesh;

template<unsigned int Dimension>
class featkAttributable {

    public:

        ~featkAttributable();

        MatrixXd getAttributeValue(size_t id) const;

    protected:

        friend class featkMesh<Dimension>;

        featkAttributable();

        void setAttribute(size_t id, std::shared_ptr<MatrixXd> attribute);

        std::map<size_t, std::shared_ptr<MatrixXd>> attributes;
};

template<unsigned int Dimension>
featkAttributable<Dimension>::featkAttributable() {

}

template<unsigned int Dimension>
featkAttributable<Dimension>::~featkAttributable() {

}


template<unsigned int Dimension>
MatrixXd featkAttributable<Dimension>::getAttributeValue(size_t id) const {

    return this->attributes.count(id) ? *(this->attributes.at(id)) : MatrixXd::Zero(1, 1);
}

template<unsigned int Dimension>
void featkAttributable<Dimension>::setAttribute(size_t id, std::shared_ptr<MatrixXd> attribute) {

    this->attributes[id] = attribute;
}


#endif // FEATKATTRIBUTUTABLE_H
