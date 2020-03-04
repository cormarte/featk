/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkBoundaryConditions.h

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
 * @class featkBoundaryConditions
 *
 * @brief Boundary conditions for finite element problems.
 *
 * featkBoundaryConditions stores the ids of the fixed degrees of freedom
 * of finite element problems along with their prescribed values.
 *
 * For efficiency reasons during the assembly process (see featkSolverBase),
 * distinction is made between zero prescrided values, for which the global
 * system vector can be trivially modified, and non-zero prescrided values,
 * involving iteration over non-zeros entries of the global system matrix.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 * @tparam Order The order of the variable the system is solved for.
 *
 */

#ifndef FEATKBOUNDARYCONDITIONS_H
#define FEATKBOUNDARYCONDITIONS_H

#include <featk/core/featkUtils.h>

#include <map>
#include <set>

template<unsigned int Dimension, unsigned int Order>
class featkBoundaryConditions {

    public:

        featkBoundaryConditions();
        ~featkBoundaryConditions();

        void addDOFValue(size_t nodeID, unsigned int nodeDOF, double value);
        std::set<size_t> getAllDOFs();
        std::map<size_t, double> getAllDOFValues();
        std::set<size_t> getNonZeroDOFs();
        std::map<size_t, double> getNonZeroDOFValues();
        std::set<size_t> getZeroDOFs();
        void setDOFValue(size_t nodeID, unsigned int nodeDOF, double value);

    private:

        std::map<size_t, double> values;
};

template<unsigned int Dimension, unsigned int Order>
featkBoundaryConditions<Dimension, Order>::featkBoundaryConditions() {

}

template<unsigned int Dimension, unsigned int Order>
featkBoundaryConditions<Dimension, Order>::~featkBoundaryConditions() {

}

template<unsigned int Dimension, unsigned int Order>
void featkBoundaryConditions<Dimension, Order>::addDOFValue(size_t nodeID, unsigned int nodeDOF, double value) {

    this->values[DOF_ID(nodeID, nodeDOF)] += value;
}

template<unsigned int Dimension, unsigned int Order>
std::set<size_t> featkBoundaryConditions<Dimension, Order>::getAllDOFs() {

    std::set<size_t> dofs;

    for (const auto& pair : this->values) {

        dofs.insert(dofs.end(), pair.first);
    }

    return dofs;
}

template<unsigned int Dimension, unsigned int Order>
std::map<size_t, double> featkBoundaryConditions<Dimension, Order>::getAllDOFValues() {

    return this->values;
}

template<unsigned int Dimension, unsigned int Order>
std::set<size_t> featkBoundaryConditions<Dimension, Order>::getNonZeroDOFs() {

    std::set<size_t> dofs;

    for (const auto& pair : this->values) {

        if (pair.second!=0.0) {

            dofs.insert(dofs.end(), pair.first);  // featkBoundaryConditions::values is already sorted by keys so insert with hint is more efficient
        }
    }

    return dofs;
}

template<unsigned int Dimension, unsigned int Order>
std::map<size_t, double> featkBoundaryConditions<Dimension, Order>::getNonZeroDOFValues() {

    std::map<size_t, double> values;

    for (const auto& pair : this->values) {

        if (pair.second!=0.0) {

            values.insert(values.end(), pair);  // featkBoundaryConditions::values is already sorted by keys so insert with hint is more efficient
        }
    }

    return values;
}

template<unsigned int Dimension, unsigned int Order>
std::set<size_t> featkBoundaryConditions<Dimension, Order>::getZeroDOFs() {

    std::set<size_t> dofs;

    for (const auto& pair : this->values) {

        if (pair.second==0.0) {

            dofs.insert(dofs.end(), pair.first);  // featkBoundaryConditions::values is already sorted by keys so insert with hint is more efficient
        }
    }

    return dofs;
}

template<unsigned int Dimension, unsigned int Order>
void featkBoundaryConditions<Dimension, Order>::setDOFValue(size_t nodeID, unsigned int nodeDOF, double value) {

    this->values[DOF_ID<Dimension, Order>(nodeID, nodeDOF)] = value;
}

#endif // FEATKBOUNDARYCONDITIONS_H
