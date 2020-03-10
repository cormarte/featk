/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkStaticSolverBase.h

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
 * @class featkStaticSolverBase
 *
 * @brief Base class for static finite element solvers.
 *
 * featkStaticSolverBase is a base class for static finite element solvers.
 *
 * featkStaticSolverBase solves finite element problems expressed as
 * \f$Kq=f\f$, where \f$K\f$ is the global system matrix returned by the
 * getGlobalSystemMatrix() function and \f$f\f$ is the global system vector
 * returned by the getGlobalSystemMatrix() function.
 *
 * Derived classes must reimplement the
 * featkSolverBase::getGlobalSystemMatrix(),
 * featkStaticSolverBase::getGlobalSystemVector(), and
 * featkSolverBase::postProcess() functions.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 * @tparam Order The order of the variable the system is solved for.
 *
 */

#ifndef FEATKSTATICSOLVERBASE_H
#define FEATKSTATICSOLVERBASE_H

#include <featk/solve/featkSolverBase.h>

template<unsigned int Dimension, unsigned int Order>
class featkStaticSolverBase : public featkSolverBase<Dimension, Order> {

    public:

        virtual ~featkStaticSolverBase();

        void solve();

    protected:

        featkStaticSolverBase();

        virtual SparseMatrix<double> getGlobalSystemMatrix()=0;
        virtual VectorXd getGlobalSystemVector()=0;
        virtual void postProcess(const VectorXd& solution)=0;
};

template<unsigned int Dimension, unsigned int Order>
featkStaticSolverBase<Dimension, Order>::featkStaticSolverBase() {

}

template<unsigned int Dimension, unsigned int Order>
featkStaticSolverBase<Dimension, Order>::~featkStaticSolverBase() {

}

template<unsigned int Dimension, unsigned int Order>
void featkStaticSolverBase<Dimension, Order>::solve() {

    SparseMatrix<double> k = this->getGlobalSystemMatrix();
    VectorXd f = this->getGlobalSystemVector();

    this->applyEBC(k, f);

    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
    //solver.setTolerance(1.0e-8);
     //solver.setMaxIterations(500);
    solver.compute(k);

    VectorXd q = solver.solve(f);

    this->postProcess(q);
}

#endif // FEATKSTATICSOLVERBASE_H
