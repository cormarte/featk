/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkLinearElasticitySolver.h

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
 * @class featkLinearElasticitySolver
 *
 * @brief Finite element solver for linear elasticity problems in Dimension
 * dimensions.
 *
 * featkLinearElasticitySolver is a finite element solver for linear
 * elasticity problems in Dimension dimensions.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 */

#ifndef FEATKLINEARELASTICITYSOLVER_H
#define FEATKLINEARELASTICITYSOLVER_H

#include <featk/solve/featkStaticSolverBase.h>

template<unsigned int Dimension>
class featkLinearElasticitySolver : public featkStaticSolverBase<Dimension, 1> {

    public:

        featkLinearElasticitySolver();
        ~featkLinearElasticitySolver();

    protected:

        SparseMatrix<double> getGlobalSystemMatrix();
        VectorXd getGlobalSystemVector();
        void postProcess(const VectorXd& u);

        std::string bodyForceAttributeName;  // Check if attributes are valid and assign their IDs to vars
        std::string outputAttributeName;
        std::string stiffnessAttributeName;
};

template<unsigned int Dimension>
featkLinearElasticitySolver<Dimension>::featkLinearElasticitySolver() {

    this->bodyForceAttributeName = "Body Force";
    this->outputAttributeName = "Displacements";
    this->stiffnessAttributeName = "Stiffness Tensor";
}

template<unsigned int Dimension>
featkLinearElasticitySolver<Dimension>::~featkLinearElasticitySolver() {

}

template<unsigned int Dimension>
SparseMatrix<double> featkLinearElasticitySolver<Dimension>::getGlobalSystemMatrix() {

    return this->getGlobalMatrixFromElements(&featkLinearElasticitySolver<Dimension>::getElementBtCBIntegralMatrix, {this->mesh->getElementAttributeID(this->stiffnessAttributeName, 4)});
}

template<unsigned int Dimension>
VectorXd featkLinearElasticitySolver<Dimension>::getGlobalSystemVector() {

    VectorXd fb = this->getGlobalVectorFromElements(&featkLinearElasticitySolver<Dimension>::getElementNtNQIntegralVector, {this->mesh->getNodeAttributeID(this->bodyForceAttributeName, 1)});
    //VectorXd fs = this->getGlobalVector(&featkStructuralAnalysisSolver<Dimension>::getElementSurfaceForceVector);
    VectorXd fs = this->getGlobalVectorFromNBCs();

    return fb+fs;
}

template<unsigned int Dimension>
void featkLinearElasticitySolver<Dimension>::postProcess(const VectorXd& u) {

    this->mesh->setNodeAttributeFromValues(this->outputAttributeName, 1, u);
}

#endif // FEATKLINEARELASTICITYSOLVER_H
