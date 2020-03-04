/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkGlobalSystemMatrixPruner.h

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
 * @class featkGlobalSystemMatrixPruner
 *
 * @brief Eigen::SparseMatrix pruner for EBC modification of the global
 * system matrix.
 *
 * featkGlobalSystemMatrixPruner is a Eigen::SparseMatrix pruner used to
 * remove non-diagonal non-zero entries of the global system matrix
 * corresponding to degrees of freedom fixed by essential boundary
 * conditions. See Eigen::SparseMatrix::prune() member function for more
 * info.
 *
 * @tparam ScalarType The scalar type of the Eigen::SparseMatrix.
 *
 */

#ifndef FEATKGLOBALSYSTEMMATRIXPRUNER_H
#define FEATKGLOBALSYSTEMMATRIXPRUNER_H

#include <set>

template<typename ScalarType>
class featkGlobalSystemMatrixPruner {

    public:

        featkGlobalSystemMatrixPruner(std::set<size_t> dofs);
        ~featkGlobalSystemMatrixPruner();

        bool operator() (const Index& row, const Index& col, const ScalarType& value) const;

    private:

        std::set<size_t> dofs;
};

template<typename ScalarType>
featkGlobalSystemMatrixPruner<ScalarType>::featkGlobalSystemMatrixPruner(std::set<size_t> dofs) {

    this->dofs = dofs;
}

template<typename ScalarType>
featkGlobalSystemMatrixPruner<ScalarType>::~featkGlobalSystemMatrixPruner() {

    this->dofs = dofs;
}

template<typename ScalarType>
bool featkGlobalSystemMatrixPruner<ScalarType>::operator() (const Index& row, const Index& col, const ScalarType& value) const {

    /* Keep diagonal elements and elements whose neither i or j index is a BC fixed dof */

    return (row == col || (this->dofs.count(row) == 0 && this->dofs.count(col) == 0));
}

#endif // FEATKGLOBALSYSTEMMATRIXPRUNER_H
