#ifndef FEATKREACTIONDIFFUSIONSOLVER_H
#define FEATKREACTIONDIFFUSIONSOLVER_H

#include <featk/solve/featkDynamicSolverBase.h>

#include <iomanip>

template<unsigned int Dimension>
class featkReactionDiffusionSolver final : public featkDynamicSolverBase<Dimension, 0> {

    public:

        featkReactionDiffusionSolver();
        ~featkReactionDiffusionSolver();

    protected:

        SparseMatrix<double> getGlobalSystemMatrix();
        virtual VectorXd getGlobalInitialVector();
        VectorXd getGlobalSystemVector(const VectorXd& u);
        void initialize();
        void intermediateProcess(const VectorXd& u, unsigned int iteration);
        void postProcess(const VectorXd& u);

        std::string diffusionAttributeName;  // Check if attributes are valid and assign their IDs to vars
        std::string inputAttributeName;
        std::string outputAttributeName;
        std::string reactionAttributeName;

        SparseMatrix<double> m;
        SparseMatrix<double> d;
        SparseMatrix<double> r;
};

template<unsigned int Dimension>
featkReactionDiffusionSolver<Dimension>::featkReactionDiffusionSolver() {

    this->diffusionAttributeName = "Diffusion Tensor";
    this->inputAttributeName = "Initial Cell Density";
    this->outputAttributeName = "Final Cell Density";
    this->reactionAttributeName = "Proliferation Rate";
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

    return this->mesh->getNodeAttributeValues(this->inputAttributeName, 0);
}

template<unsigned int Dimension>
VectorXd featkReactionDiffusionSolver<Dimension>::getGlobalSystemVector(const VectorXd &u) {

    /*
        The last term of the equation should be integral(Nt(NU)^2) instead of integral(NtN)*U^2 = MU^2.
        See Mocenni et al. 2011. for handling of polynomial reaction terms in FEM.
    */

    size_t id = this->mesh->setNodeAttributeFromValues("tmp", 0, u);
    VectorXd ru2 = this->getGlobalVectorFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtCNQNQIntegralVector, {this->mesh->getElementAttributeID(this->reactionAttributeName, 0), id});
    return this->m*u + this->timeStep*(this->r*u-ru2);

    //return this->m*u + this->timeStep*this->r*(u-u.cwiseProduct(u));
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::initialize() {

    this->m = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtNIntegralMatrix, {});
    cout << "featkReactionDiffusionSolver: Info: M matrix assembled." << endl;
    this->d = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementBtCBIntegralMatrix, {this->mesh->getElementAttributeID(this->diffusionAttributeName, 2)});
    cout << "featkReactionDiffusionSolver: Info: D matrix assembled." << endl;
    this->r = this->getGlobalMatrixFromElements(&featkReactionDiffusionSolver<Dimension>::getElementNtCNIntegralMatrix, {this->mesh->getElementAttributeID(this->reactionAttributeName, 0)});
    cout << "featkReactionDiffusionSolver: Info: R matrix assembled." << endl;
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::intermediateProcess(const VectorXd &u, unsigned int iteration) {

    std::ostringstream stream;
    stream << fixed << std::setprecision(2) << iteration*this->timeStep;

    this->mesh->setNodeAttributeFromValues(this->outputAttributeName + " (" + stream.str() + ")", 0, u);
}

template<unsigned int Dimension>
void featkReactionDiffusionSolver<Dimension>::postProcess(const VectorXd& u) {

    this->mesh->setNodeAttributeFromValues(this->outputAttributeName, 0, u);
    this->mesh->computeNodeBQ<0>(this->outputAttributeName, this->outputAttributeName + " Gradient");
}

#endif // FEATKREACTIONDIFFUSIONSOLVER_H
