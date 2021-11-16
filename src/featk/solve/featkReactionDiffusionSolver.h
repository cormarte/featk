/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkReactionDiffusionSolver.h

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

 * @class featkReactionDiffusionSolver
 *
 * @brief Dynamic finite element solver for reaction-diffusion problems in
 * Dimension dimensions.
 *
 * featkReactionDiffusionSolver is a dynamic finite element solver for
 * reaction-diffusion problems with a logistic reaction term in Dimension
 * dimensions.
 *
 * featkReactionDiffusionSolver solves
 *
 * \f[
 *
 * \begin{equation}
 * \begin{cases}
 * & \frac{\partial u(\bar{r}, t)}{\partial t} &= \bar{\nabla} \cdot \left( \bar{\bar{D}}(\bar{r}) \bar{\nabla} u(\bar{r}, t) \right) + \rho(\bar{r}) u(\bar{r}, t) \left( 1-u(\bar{r}, t) \right) \\
 * & u(\bar{r}, 0) &= u_0(\bar{r}) \quad \forall \bar{r} \in \Omega \\
 * & \left( \bar{\bar{D}}(\bar{r}) \bar{\nabla} u(\bar{r}, t) \right) \cdot \bar{n} &= 0 \quad \forall \bar{r} \in \partial \Omega
 * \end{cases}
 * \end{equation}
 *
 * \f]
 *
 * where \f$u(\bar{r}, t)\f$ is a normalized function of cartesian position
 * \f$\bar{r}\f$ and time \f$t\f$ such that \f$u: \Omega \times \left[ 0, T
 * \right] \to \left[ 0, 1 \right]; (\bar{r}, t) \mapsto u(\bar{r}, t)\f$
 * where \f$\Omega \subset \mathbb{R}^{N}\f$ is the solving spatial domain
 * and \f$T\f$ is the simulated time; \f$\bar{\bar{D}}(\bar{r})\f$ and
 * \f$\rho(\bar{r})\f$ are respectively the value of the diffusion tensor
 *  and proliferation rate at position \f$\bar{r}\f$; \f$\partial\Omega\f$
 * is the boundary of \f$\Omega\f$; and \f$\bar{n}\f$ is the unit normal
 * vector on \f$\partial \Omega\f$.
 *
 * @tparam Dimension The cartesian dimension of the problem.
 *
 */

#ifndef FEATKREACTIONDIFFUSIONSOLVER_H
#define FEATKREACTIONDIFFUSIONSOLVER_H

#include <featk/solve/featkDynamicSolverBase.h>

#include <iomanip>

template<unsigned int Dimension>
class featkReactionDiffusionSolver final : public featkDynamicSolverBase<Dimension, 0> {

    public:

        featkReactionDiffusionSolver();
        ~featkReactionDiffusionSolver();

        void setDiffusionElementAttributeName(std::string name);  // Check if attributes are valid and assign their IDs to vars
        void setInputNodeAttributeName(std::string name);
        void setOutputNodeAttributeName(std::string name);
        void setReactionElementAttributeName(std::string name);
        void setUseSpeedHack(bool use);

    protected:

        SparseMatrix<double> getGlobalSystemMatrix();
        VectorXd getGlobalInitialVector();
        VectorXd getGlobalSystemVector(const VectorXd& u);
        void initialize();
        void intermediateProcess(const VectorXd& u, unsigned int iteration);
        void postProcess(const VectorXd& u);

        std::string diffusionElementAttributeName;  // Check if attributes are valid and assign their IDs to vars
        std::string inputNodeAttributeName;
        std::string outputNodeAttributeName;
        std::string reactionElementAttributeName;
        bool useSpeedHack;

        SparseMatrix<double> m;
        SparseMatrix<double> d;
        SparseMatrix<double> r;
};

template<unsigned int Dimension>
featkReactionDiffusionSolver<Dimension>::featkReactionDiffusionSolver() {

    this->diffusionElementAttributeName = "Diffusion Tensor";
    this->inputNodeAttributeName = "Initial Cell Density";
    this->outputNodeAttributeName = "Final Cell Density";
    this->reactionElementAttributeName = "Proliferation Rate";
    this->useSpeedHack = true;
}

template<unsigned int Dimension>
featkReactionDiffusionSolver<Dimension>::~featkReactionDiffusionSolver() {

}

template<unsigned int Dimension>
SparseMatrix<double> featkReactionDiffusionSolver<Dimension>::getGlobalSystemMatrix() {

    return this->m + this->timeStep*d;
}

template<unsigned int Dimension>
VectorXd featkReactionDiffusionSolver<Dimension>::getGlobalInitialVector() {

    return this->mesh->getNodeAttributeValues(this->inputNodeAttributeName, 0);
}

template<unsigned int Dimension>
VectorXd featkReactionDiffusionSolver<Dimension>::getGlobalSystemVector(const VectorXd &u) {

    /*
        The last term of the equation should be integral(Nt(NU)^2) instead of integral(NtN)*U^2 = MU^2.
        See Mocenni et al. 2011. for handling of polynomial reaction terms in FEM.
    */

    VectorXd f;

    if (this->useSpeedHack) {

        f = this->m*u + this->timeStep*this->r*(u-u.cwiseProduct(u));
    }

    else {

        size_t id = this->mesh->setNodeAttributeFromValues("tmp", 0, u);
        VectorXd ru2 = this->getGlobalVectorFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtCNQNQIntegralVector, {this->mesh->getElementAttributeID(this->reactionElementAttributeName, 0), id});
        f = this->m*u + this->timeStep*(this->r*u-ru2);
    }

    return f;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::initialize() {

    this->m = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtNIntegralMatrix, {});
    cout << "featkReactionDiffusionSolver: Info: M matrix assembled." << endl;
    this->d = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementBtCBIntegralMatrix, {this->mesh->getElementAttributeID(this->diffusionElementAttributeName, 2)});
    cout << "featkReactionDiffusionSolver: Info: D matrix assembled." << endl;
    this->r = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtCNIntegralMatrix, {this->mesh->getElementAttributeID(this->reactionElementAttributeName, 0)});
    cout << "featkReactionDiffusionSolver: Info: R matrix assembled." << endl;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::intermediateProcess(const VectorXd &u, unsigned int iteration) {

    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << iteration*this->timeStep;

    this->mesh->setNodeAttributeFromValues(this->outputNodeAttributeName + " (" + stream.str() + ")", 0, u);
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::postProcess(const VectorXd& u) {

    this->mesh->setNodeAttributeFromValues(this->outputNodeAttributeName, 0, u);
    // this->mesh->computeNodeBQ<0>(this->outputNodeAttributeName, this->outputNodeAttributeName + " Gradient");  // No more perfmored here since gradient is zero along CSF boundaries
}


template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::setDiffusionElementAttributeName(std::string name) {

    this->diffusionElementAttributeName = name;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::setInputNodeAttributeName(std::string name) {

    this->inputNodeAttributeName = name;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::setOutputNodeAttributeName(std::string name) {

    this->outputNodeAttributeName = name;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::setReactionElementAttributeName(std::string name) {

    this->reactionElementAttributeName = name;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::setUseSpeedHack(bool use) {

    this->useSpeedHack = use;
}

#endif // FEATKREACTIONDIFFUSIONSOLVER_H
