/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkSolverBase.h

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

 * @class featkSolverBase
 *
 * @brief Base class for finite element solvers in Dimension dimensions.
 *
 * featkReactionDiffusionSolver is base class for finite element solvers in
 * Dimension dimensions. featkReactionDiffusionSolver gathers featkMesh and
 * featkBoundaryCondition objects involved in a finite element problem.
 *
 * featkReactionDiffusionSolver provides functions for assembling global
 * matrices and vectors from featkElementInterface quantities getters as
 * well as for applyingy essential boundary conditions to the global system
 * matrix and vector.
 *
 * featkReactionDiffusionSolver also defines the update() routine as the
 * sucession of calls to initialize() and solve() functions.
 * The initialize() function is used to initialize matrices and vectors that
 * are not subject to change throughout the solving process. The solve()
 * function solves the probleù itself, at once or iteratively. The
 * postProcess() function is called at the end of the solve routine to
 * assign solution and derivative quantities to the input featkMesh.
 *
 * Derived classes must reimplement the featkSolverBase::solve(),
 * featkSolverBase::getGlobalSystemMatrix(), and
 * featkSolverBase::postProcess() functions and may reimplement the
 * featkSolverBase::initialize() function if needed.
 *
 * @todo Make featkSolverBase inherit from featkMeshConsumerBase.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 * @tparam Order The order of the variable the system is solved for.
 *
 */

#ifndef FEATKSOLVERBASE_H
#define FEATKSOLVERBASE_H

#include <featk/core/featkUtils.h>
#include <featk/geometry/featkMesh.h>
#include <featk/solve/featkBoundaryConditions.h>
#include <featk/solve/featkGlobalSystemMatrixPruner.h>

#include <Eigen/Sparse>
#include <iostream>

using namespace Eigen;

template<unsigned int Dimension, unsigned int Order>
class featkSolverBase {

    public:

        virtual ~featkSolverBase();

        virtual void initialize();
        virtual void solve()=0;

        void update();
        void setEssentialBoundaryConditions(featkBoundaryConditions<Dimension, Order>* conditions);
        void setInputMesh(featkMesh<Dimension>* mesh);
        void setNaturalBoundaryConditions(featkBoundaryConditions<Dimension, Order>* conditions);

        static const unsigned int dofsPerNode = POWER(Dimension, Order);

    protected:

        static MatrixXd getElementBtBIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static MatrixXd getElementBtCBIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static MatrixXd getElementNtCNIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static VectorXd getElementNtCNQNQIntegralVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static MatrixXd getElementNtNIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static VectorXd getElementNtNQIntegralVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);
        static VectorXd getElementQVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs);

        featkSolverBase();

        virtual SparseMatrix<double> getGlobalSystemMatrix()=0;
        virtual void postProcess(const VectorXd& solution)=0;

        void applyEBC(SparseMatrix<double>& k, VectorXd& f);
        void applyEBCToGlobalSystemVector(const SparseMatrix<double>& globalStiffnessMatrix, VectorXd& f);
        void applyEBCToGlobalSystemMatrix(SparseMatrix<double>& k);
        VectorXd getEBCModifiedGlobalSystemVector(const SparseMatrix<double>& k, const VectorXd& vector);
        SparseMatrix<double> getEBCModifiedGlobalSystemMatrix(const SparseMatrix<double>& matrix);
        SparseMatrix<double> getGlobalMatrixFromElements(MatrixXd (*getElementMatrix)(featkElementInterface<Dimension>*, std::vector<size_t>), std::vector<size_t> attributeIDs);           // Assembles global matrix from element matrix getter
        // void getGlobalMatrixFromElements(MatrixXd (*getElementMatrix)(featkElementInterface<Dimensions>*), SparseMatrix<double>& k);  // Check if performs faster (i.e. if NRVO is not applied to Eigen::SparseMatrix)
        VectorXd getGlobalVectorFromNBCs();
        VectorXd getGlobalVectorFromElements(VectorXd (*getElementVector)(featkElementInterface<Dimension>*, std::vector<size_t>), std::vector<size_t> attributeIDs);                        // Assembles global vector from element vector getter

        featkBoundaryConditions<Dimension, Order>* essentialBoundaryConditions;
        featkMesh<Dimension>* mesh;
        featkBoundaryConditions<Dimension, Order>* naturalBoundaryConditions;
        size_t numberOfDOFs;
};

template<unsigned int Dimension, unsigned int Order>
MatrixXd featkSolverBase<Dimension, Order>::getElementBtBIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getBtBIntegralMatrix<Order>();
}

template<unsigned int Dimension, unsigned int Order>
MatrixXd featkSolverBase<Dimension, Order>::getElementBtCBIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getBtCBIntegralMatrix<Order>(attributeIDs[0]);
}

template<unsigned int Dimension, unsigned int Order>
MatrixXd featkSolverBase<Dimension, Order>::getElementNtCNIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getNtCNIntegralMatrix<Order>(attributeIDs[0]);
}

template<unsigned int Dimension, unsigned int Order>
MatrixXd featkSolverBase<Dimension, Order>::getElementNtNIntegralMatrix(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getNtNIntegralMatrix<Order>();
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getElementNtNQIntegralVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getNtNQIntegralVector<Order>(attributeIDs[0]);
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getElementNtCNQNQIntegralVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getNtCNQNQIntegralVector<Order>(attributeIDs[0], attributeIDs[1]);
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getElementQVector(featkElementInterface<Dimension>* element, std::vector<size_t> attributeIDs) {

    return element->getQVector<Order>(attributeIDs[0]);
}


template<unsigned int Dimension, unsigned int Order>
featkSolverBase<Dimension, Order>::featkSolverBase() {

    this-> essentialBoundaryConditions = nullptr;
    this->mesh = nullptr;
    this->naturalBoundaryConditions = nullptr;
    this->numberOfDOFs = 0;
}

template<unsigned int Dimension, unsigned int Order>
featkSolverBase<Dimension, Order>::~featkSolverBase() {

}


template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::applyEBC(SparseMatrix<double>& k, VectorXd& f) {

    this->applyEBCToGlobalSystemVector(k, f);
    this->applyEBCToGlobalSystemMatrix(k);
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::applyEBCToGlobalSystemVector(const SparseMatrix<double>& k, VectorXd& f) {

    std::map<size_t, double> allDOFValues = this->essentialBoundaryConditions->getAllDOFValues();
    std::map<size_t, double> nonZeroDOFValues = this->essentialBoundaryConditions->getNonZeroDOFValues();
    std::cout << "featkSolverBase: Info: System has " << nonZeroDOFValues.size() << " non-zero essential boundary conditions." << std::endl;

    for (const auto& pair : allDOFValues) {

        f(pair.first, 0) = pair.second;
    }

    if (nonZeroDOFValues.size()!=0) {

        for (size_t i=0; i!=k.outerSize(); i++) {  // Iterate over all non zero entries of the sparse matrix

            for (SparseMatrix<double>::InnerIterator it(k, i); it; ++it) {

                if (allDOFValues.count(it.row()) == 0 && nonZeroDOFValues.count(it.col()) != 0) {

                    f(it.row(), 0) -= it.value()*nonZeroDOFValues[it.col()];
                }
            }
        }
    }
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::applyEBCToGlobalSystemMatrix(SparseMatrix<double>& k) {

    std::set<size_t> allDOFs = this->essentialBoundaryConditions->getAllDOFs();
    cout << "featkSolverBase: Info: System has " << allDOFs.size() << " essential boundary conditions." << endl;

    for (int i=0; i!=k.outerSize(); i++) {

        for (SparseMatrix<double>::InnerIterator it(k, i); it; ++it) {

            if (allDOFs.count(it.row()) != 0 && it.row() == it.col()) {  // Diagonal elements are always non zero or explicit zero, so first test if dof is fixed.

                it.valueRef() = 1.0;
            }
        }
    }

    /* The global stiffness matrix diagonal elements are non-zero or at least an explicit zero, so it is better
     * to first change their value to 1, and prune non-diagonal elements afterwards. Indeed, non-zero element
     * insertion into an Eigen::SparseMatrix is very expensive */

    k.prune(featkGlobalSystemMatrixPruner<double>(allDOFs));
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getEBCModifiedGlobalSystemVector(const SparseMatrix<double>& k, const VectorXd& vector) {

    VectorXd f = vector;
    this->applyEBCToGlobalSystemVector(k, f);

    return f;
}

template<unsigned int Dimension, unsigned int Order>
SparseMatrix<double> featkSolverBase<Dimension, Order>::getEBCModifiedGlobalSystemMatrix(const SparseMatrix<double>& matrix) {

    SparseMatrix<double> k = matrix;
    this->applyEBCToGlobalSystemMatrix(k);

    return k;
}

template<unsigned int Dimension, unsigned int Order>
SparseMatrix<double> featkSolverBase<Dimension, Order>::getGlobalMatrixFromElements(MatrixXd (*getElementMatrix)(featkElementInterface<Dimension>*, std::vector<size_t>), std::vector<size_t> attributeIDs) {

    std::map<size_t, std::map<size_t, double>> coefficients;

    for (featkElementInterface<Dimension>* element : this->mesh->getElements()) {

        std::vector<size_t> elementDOFs;

        for (featkNode<Dimension>* node : element->getNodes()) {

            for (unsigned int dof=0; dof!=this->dofsPerNode; dof++) {

                elementDOFs.push_back(DOF_ID<Dimension, Order>(node->getID(), dof));
            }
        }

        MatrixXd elementMatrix = getElementMatrix(element, attributeIDs);

        for (int i=0; i!=elementDOFs.size(); i++) {

            for (int j=0; j!=elementDOFs.size(); j++) {

                size_t k = elementDOFs[i];
                size_t l = elementDOFs[j];

                coefficients[k][l] += elementMatrix(i, j);
            }
        }
    }

    std::vector<Triplet<double, size_t>> triplets;

    for (auto const &a : coefficients) {

        size_t i = a.first;

        for (auto const &b : a.second) {

            size_t j = b.first;
            triplets.push_back(Triplet<double, size_t>(i, j, b.second));
        }
    }

    SparseMatrix<double> k(this->numberOfDOFs, this->numberOfDOFs);
    k.setFromTriplets(triplets.begin(), triplets.end());

    return k;  // Make sure NRVO is applied here to avoid copying a huge Eigen::SparseMatrix
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getGlobalVectorFromNBCs() {

    VectorXd f = VectorXd::Zero(this->numberOfDOFs);

    for (const auto& pair : this->naturalBoundaryConditions->getNonZeroDOFValues()) {

        f(pair.first, 0) = pair.second;
    }

    return f;
}

template<unsigned int Dimension, unsigned int Order>
VectorXd featkSolverBase<Dimension, Order>::getGlobalVectorFromElements(VectorXd (*getElementVector)(featkElementInterface<Dimension>*, std::vector<std::size_t>), std::vector<size_t> attributeIDs) {

    VectorXd f = VectorXd::Zero(this->numberOfDOFs);

    for (featkElementInterface<Dimension>* element : this->mesh->getElements()) {

        VectorXd elementVector = getElementVector(element, attributeIDs);
        size_t i = 0;

        for (featkNode<Dimension>* node : element->getNodes()) {

            for (unsigned int dof=0; dof!=this->dofsPerNode; dof++) {

                f(DOF_ID<Dimension, Order>(node->getID(), dof), 0) += elementVector(i, 0);
                i++;
            }
        }
    }

    return f;
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::initialize() {

}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::update() {

    this->initialize();
    this->solve();
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::setEssentialBoundaryConditions(featkBoundaryConditions<Dimension, Order>* conditions) {

    this->essentialBoundaryConditions = conditions;
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::setInputMesh(featkMesh<Dimension>* mesh) {

    this->mesh = mesh;
    this->numberOfDOFs = mesh->getNumberOfNodes()*this->dofsPerNode;
}

template<unsigned int Dimension, unsigned int Order>
void featkSolverBase<Dimension, Order>::setNaturalBoundaryConditions(featkBoundaryConditions<Dimension, Order>* conditions) {

    this->naturalBoundaryConditions = conditions;
}

#endif // FEATKSOLVERBASE_H
