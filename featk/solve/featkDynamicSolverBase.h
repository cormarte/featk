/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkDynamicSolverBase.h

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
 * @class featkDynamicSolverBase
 *
 * @brief Base class for dynamic finite element solvers.
 *
 * featkDynamicSolverBase is a base class for dynamic finite element solvers.
 * featkDynamicSolverBase introduces the concept of iterative time-step
 * solving and provides convenient mechanisms for stepwise global system
 * vector recalculation and intermediate solution processing.
 *
 * Typical iterative solving is perfomed as follows:
 *
 * 1) Global matrices involved that are not subject to change throughout the
 * solving process are computed inside the overriden virtual member function
 * featkSolverBase::initialize().
 *
 * 2) The global system matrix is generated from the global matrices
 * assembled in 1) and EBCs are applied to it.
 *
 * 3) An initial solution vector is generated.
 *
 * 4) For each iteration:\n
 * a. A global system vector is generated from the global matrices computed in 1) and the previous solution.\n
 * b. EBCs are applied to the global vector.\n
 * c. System is solved.\n
 * d. Current solution is cut off if required.\n
 * e. Intermediate processing is performed on the current solution if required.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 * @tparam Order The order of the variable the system is solved for.
 *
 */

#ifndef FEATKDYNAMICSOLVERBASE_H
#define FEATKDYNAMICSOLVERBASE_H

#include <featk/solve/featkSolverBase.h>

template<unsigned int Dimension, unsigned int Order>
class featkDynamicSolverBase : public featkSolverBase<Dimension, Order> {

    public:

        virtual ~featkDynamicSolverBase();

        void solve();

        void setDoCutoff(bool doCutoff);
        void setDoLowerCutoff(bool doCutoff);
        void setDoUpperCutoff(bool doCutoff);
        void setIntermediateProcessIterations(std::vector<unsigned int> iterations);
        void setLowerCutoffValue(double value);
        void setNumberOfIterations(unsigned int iterations);
        void setTimeStep(double step);
        void setUpperCutoffValue(double value);

    protected:

        featkDynamicSolverBase();

        virtual SparseMatrix<double> getGlobalSystemMatrix()=0;
        virtual VectorXd getGlobalInitialVector()=0;
        virtual VectorXd getGlobalSystemVector(const VectorXd& u)=0;
        virtual void intermediateProcess(const VectorXd& u, unsigned int iteration);
        virtual void postProcess(const VectorXd& solution)=0;

        bool doLowerCutoff;
        bool doUpperCutoff;
        double lowerCutoffValue;
        double upperCutoffValue;

        unsigned int numberOfIterations;
        double timeStep;
        std::vector<unsigned int> intermediateProcessIterations;
};

template<unsigned int Dimension, unsigned int Order>
featkDynamicSolverBase<Dimension, Order>::featkDynamicSolverBase() {

    this->doLowerCutoff = false;
    this->doUpperCutoff = false;
    this->lowerCutoffValue = 0.0;
    this->upperCutoffValue = 0.0;

    this->numberOfIterations = 500;
    this->timeStep = 1.0;
}

template<unsigned int Dimension, unsigned int Order>
featkDynamicSolverBase<Dimension, Order>::~featkDynamicSolverBase() {

}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::intermediateProcess(const VectorXd &u, unsigned int iteration) {

}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setDoCutoff(bool doCutoff) {

    this->doLowerCutoff = doCutoff;
    this->doUpperCutoff = doCutoff;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setDoLowerCutoff(bool doCutoff) {

    this->doLowerCutoff = doCutoff;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setDoUpperCutoff(bool doCutoff) {

    this->doUpperCutoff = doCutoff;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setIntermediateProcessIterations(std::vector<unsigned int> iterations) {

    this->intermediateProcessIterations = iterations;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setLowerCutoffValue(double value) {

    this->lowerCutoffValue = value;
    this->doLowerCutoff = true;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setNumberOfIterations(unsigned int iterations) {

    this->numberOfIterations = iterations;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setTimeStep(double step) {

    this->timeStep = step;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::setUpperCutoffValue(double value) {

    this->upperCutoffValue	= value;
    this->doUpperCutoff = true;
}

template<unsigned int Dimension, unsigned int Order>
void featkDynamicSolverBase<Dimension, Order>::solve() {

    cout << "featkDynamicSolverBase: Info: System has " << this->numberOfDOFs << " degrees of freedom." << endl;

    SparseMatrix<double> globalSystemMatrix = this->getGlobalSystemMatrix();
    SparseMatrix<double> k = this->getEBCModifiedGlobalSystemMatrix(globalSystemMatrix);

    cout << "featkDynamicSolverBase: Info: System matrix density is " << globalSystemMatrix.nonZeros() << "/" << this->numberOfDOFs*this->numberOfDOFs << "." << endl;

    /*k.makeCompressed();

    SimplicialLLT<SparseMatrix<double>> solver;
    solver.compute(k);

    if (solver.info() != Success) {

        cout << "featkDynamicSolverBase: Warning: Global system matrix decomposition failed." << endl;
    }*/

    //BiCGSTAB<SparseMatrix<double, RowMajor>> solver;  // OpenMP parallelized only for RowMajor. BiCGSTAB is more general than CG and works for all kind of matrices.
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;  // Only for symmetric positive definite matrices, a bit faster than BiCGSTAB in this case.
    solver.compute(k);

    VectorXd u = this->getGlobalInitialVector();
    VectorXd f = VectorXd(this->numberOfDOFs);

    for (unsigned int i=0; i!=this->numberOfIterations; i++) {

        f = this->getGlobalSystemVector(u);
        this->applyEBCToGlobalSystemVector(globalSystemMatrix, f);
        u = solver.solveWithGuess(f, u).head(this->numberOfDOFs);
        //u = solver.solve(u).head(this->numberOfDOFs);

        if (this->doLowerCutoff) {

            u = (u.array() < this->lowerCutoffValue).select(this->lowerCutoffValue, u);
        }

        if (this->doUpperCutoff) {

            u = (u.array() > this->upperCutoffValue).select(this->upperCutoffValue, u);
        }

        cout << "featkDynamicSolverBase: Info: Iteration " << i+1 << "/" << this->numberOfIterations << " solved (" << solver.iterations() << " iterations, error: " << solver.error() << ")." << endl;
        //cout << "featkDynamicSolverBase: Info: Iteration " << i+1 << "/" << this->numberOfIterations << " solved." << endl;

        if (find(this->intermediateProcessIterations.begin(), this->intermediateProcessIterations.end(), i) != this->intermediateProcessIterations.end()) {

            this->intermediateProcess(u, i);
        }
    }

    cout << "featkDynamicSolverBase: Info: System solved" << endl;

    this->postProcess(u);

    cout << "featkDynamicSolverBase: Info: Post processing done" << endl;
}

#endif // FEATKDYNAMICSOLVERBASE_H
