#ifndef FEATK2POPULATIONSREACTIONDIFFUSIONSOLVER_H
#define FEATK2POPULATIONSREACTIONDIFFUSIONSOLVER_H

#include <featk/solve/featkDynamicSolverBase.h>

#include <iomanip>

template<unsigned int Dimension>
class featk2PopulationsReactionDiffusionSolver final : public featkDynamicSolverBase<Dimension, 0> {

    public:

        featk2PopulationsReactionDiffusionSolver();
        ~featk2PopulationsReactionDiffusionSolver();

        void solve();

        void setUseSpeedHack(bool use);

    protected:

        SparseMatrix<double> getGlobalSystemMatrix();
        VectorXd getGlobalInitialVector();
        VectorXd getGlobalSystemVector(const VectorXd& u);
        void initialize();
        void intermediateProcess(const VectorXd& u, unsigned int iteration);
        void postProcess(const VectorXd& u);


        SparseMatrix<double> getGlobalSystemMatrix(double df);
        VectorXd getGlobalInitialVector(std::string inputNodeAttributeName);
        VectorXd getGlobalSystemVector(const VectorXd& u, const VectorXd& t, double rf);
        VectorXd getGlobalSystemVector(const VectorXd& u, const VectorXd& t, double df, double pf);
        void intermediateProcess(const VectorXd& u, std::string name, unsigned int iteration);
        void postProcess(const VectorXd& u, std::string name);

        std::string diffusionElementAttributeName;  // Check if attributes are valid and assign their IDs to vars
        std::string inputNodeAttributeName1;
        std::string inputNodeAttributeName2;
        std::string outputNodeAttributeName1;
        std::string outputNodeAttributeName2;
        std::string reactionElementAttributeName;
        bool useSpeedHack;

        SparseMatrix<double> m;
        SparseMatrix<double> d;
        SparseMatrix<double> r;
};

template<unsigned int Dimension>
featk2PopulationsReactionDiffusionSolver<Dimension>::featk2PopulationsReactionDiffusionSolver() {

    this->diffusionElementAttributeName = "Diffusion Tensor";
    this->inputNodeAttributeName1 = "Initial Cell Density 1";
    this->inputNodeAttributeName2 = "Initial Cell Density 2";
    this->outputNodeAttributeName1 = "Final Cell Density 1";
    this->outputNodeAttributeName2 = "Final Cell Density 2";
    this->reactionElementAttributeName = "Proliferation Rate";
    this->useSpeedHack = true;
}

template<unsigned int Dimension>
featk2PopulationsReactionDiffusionSolver<Dimension>::~featk2PopulationsReactionDiffusionSolver() {

}

template<unsigned int Dimension>
SparseMatrix<double> featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalSystemMatrix() {

    return this->m;
}

template<unsigned int Dimension>
VectorXd featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalInitialVector() {

    return VectorXd(1, 1);
}

template<unsigned int Dimension>
VectorXd featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalSystemVector(const VectorXd &u) {

    return VectorXd(1, 1);
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::initialize() {

    this->m = this->getGlobalMatrixFromElements(&featk2PopulationsReactionDiffusionSolver<Dimension>::getElementNtNIntegralMatrix, {});
    cout << "featk2PopulationsReactionDiffusionSolver: Info: M matrix assembled." << endl;
    this->d = this->getGlobalMatrixFromElements(&featk2PopulationsReactionDiffusionSolver<Dimension>::getElementBtCBIntegralMatrix, {this->mesh->getElementAttributeID(this->diffusionElementAttributeName, 2)});
    cout << "featk2PopulationsReactionDiffusionSolver: Info: D matrix assembled." << endl;
    this->r = this->getGlobalMatrixFromElements(&featk2PopulationsReactionDiffusionSolver<Dimension>::getElementNtCNIntegralMatrix, {this->mesh->getElementAttributeID(this->reactionElementAttributeName, 0)});
    cout << "featk2PopulationsReactionDiffusionSolver: Info: R matrix assembled." << endl;
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::intermediateProcess(const VectorXd &u, unsigned int iteration) {

}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::postProcess(const VectorXd& u) {

}


template<unsigned int Dimension>
VectorXd featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalInitialVector(std::string name) {

    return this->mesh->getNodeAttributeValues(name, 0);
}

template<unsigned int Dimension>
SparseMatrix<double> featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalSystemMatrix(double df) {

    return this->m + df*this->timeStep*d;
}

template<unsigned int Dimension>
VectorXd featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalSystemVector(const VectorXd &u, const VectorXd &t, double rf) {

    /*
        The last term of the equation should be integral(Nt(NU)^2) instead of integral(NtN)*U^2 = MU^2.
        See Mocenni et al. 2011. for handling of polynomial reaction terms in FEM.
    */

    VectorXd f;

    if (this->useSpeedHack) {

        f = this->m*u + rf*this->timeStep*this->r*(u-u.cwiseProduct(t));
    }

    /*else {

        size_t id = this->mesh->setNodeAttributeFromValues("tmp", 0, u);
        VectorXd ru2 = this->getGlobalVectorFromElements(&featk2PopulationsReactionDiffusionSolver<Dimension>::getElementNtCNQNQIntegralVector, {this->mesh->getElementAttributeID(this->reactionElementAttributeName, 0), id});
        f = this->m*u + this->timeStep*(this->r*u-ru2);
    }*/

    return f;
}


template<unsigned int Dimension>
VectorXd featk2PopulationsReactionDiffusionSolver<Dimension>::getGlobalSystemVector(const VectorXd &u, const VectorXd &t, double df, double pf) {

    /*
        The last term of the equation should be integral(Nt(NU)^2) instead of integral(NtN)*U^2 = MU^2.
        See Mocenni et al. 2011. for handling of polynomial reaction terms in FEM.
    */

    VectorXd f;

    if (this->useSpeedHack) {

        f = this->m*u -df*this->timeStep*this->d*t + pf*this->timeStep*this->r*(u-u.cwiseProduct(t));
    }

    /*else {

        size_t id = this->mesh->setNodeAttributeFromValues("tmp", 0, u);
        VectorXd ru2 = this->getGlobalVectorFromElements(&featk2PopulationsReactionDiffusionSolver<Dimension>::getElementNtCNQNQIntegralVector, {this->mesh->getElementAttributeID(this->reactionElementAttributeName, 0), id});
        f = this->m*u + this->timeStep*(this->r*u-ru2);
    }*/

    return f;
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::intermediateProcess(const VectorXd &u, std::string name, unsigned int iteration) {

    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << iteration*this->timeStep;

    this->mesh->setNodeAttributeFromValues(name + " (" + stream.str() + ")", 0, u);
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::setUseSpeedHack(bool use) {

    this->useSpeedHack = use;
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::postProcess(const VectorXd& u, std::string name) {

    this->mesh->setNodeAttributeFromValues(name, 0, u);
    this->mesh->computeNodeBQ<0>(name, name + " Gradient");
}

template<unsigned int Dimension>
void featk2PopulationsReactionDiffusionSolver<Dimension>::solve() {

    cout << "featk2PopulationsReactionDiffusionSolver: Info: System has " << this->numberOfDOFs << " degrees of freedom." << endl;

    //SparseMatrix<double> globalSystemMatrix1 = this->getGlobalSystemMatrix(1.0);
    //SparseMatrix<double> globalSystemMatrix2 = this->getGlobalSystemMatrix(4.0);  // Population 2 diffuses 2 times faster
    SparseMatrix<double> globalSystemMatrix1 = this->m;
    SparseMatrix<double> globalSystemMatrix2 = this->m;
    SparseMatrix<double> k1 = this->getEBCModifiedGlobalSystemMatrix(globalSystemMatrix1);
    SparseMatrix<double> k2 = this->getEBCModifiedGlobalSystemMatrix(globalSystemMatrix2);

    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver1;  // Only for symmetric positive definite matrices, a bit faster than BiCGSTAB in this case.
    solver1.compute(k1);

    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver2;  // Only for symmetric positive definite matrices, a bit faster than BiCGSTAB in this case.
    solver2.compute(k2);

    VectorXd u1 = this->getGlobalInitialVector(this->inputNodeAttributeName1);
    VectorXd u2 = this->getGlobalInitialVector(this->inputNodeAttributeName2);

    VectorXd f1 = VectorXd(this->numberOfDOFs);
    VectorXd f2 = VectorXd(this->numberOfDOFs);

    for (unsigned int i=0; i!=this->numberOfIterations; i++) {

        VectorXd t = u1 + u2;

        /*f1 = this->getGlobalSystemVector(u1, t, 1.5);
        f2 = this->getGlobalSystemVector(u2, t, 1.0);*/

        f1 = this->getGlobalSystemVector(u1, t, 1.0, 1.5);
        f2 = this->getGlobalSystemVector(u2, t, 2.0, 1.0);

        this->applyEBCToGlobalSystemVector(globalSystemMatrix1, f1);
        this->applyEBCToGlobalSystemVector(globalSystemMatrix2, f2);

        u1 = solver1.solveWithGuess(f1, u1).head(this->numberOfDOFs);
        cout << "featk2PopulationsReactionDiffusionSolver: Info: Population 1, Iteration " << i+1 << "/" << this->numberOfIterations << " solved (" << solver1.iterations() << " iterations, error: " << solver1.error() << ")." << endl;

        u2 = solver2.solveWithGuess(f2, u2).head(this->numberOfDOFs);
        cout << "featk2PopulationsReactionDiffusionSolver: Info: Population 2, Iteration " << i+1 << "/" << this->numberOfIterations << " solved (" << solver2.iterations() << " iterations, error: " << solver2.error() << ")." << endl;

        if (this->doLowerCutoff) {

            u1 = (u1.array() < this->lowerCutoffValue).select(this->lowerCutoffValue, u1);
            u2 = (u2.array() < this->lowerCutoffValue).select(this->lowerCutoffValue, u2);
        }

        if (this->doUpperCutoff) {

            u1 = (u1.array() > this->upperCutoffValue).select(this->upperCutoffValue, u1);
            u2 = (u2.array() > this->upperCutoffValue).select(this->upperCutoffValue, u2);
        }

        if (find(this->intermediateProcessIterations.begin(), this->intermediateProcessIterations.end(), i) != this->intermediateProcessIterations.end()) {

            this->intermediateProcess(u1, i);
        }
    }

    cout << "featk2PopulationsReactionDiffusionSolver: Info: System solved" << endl;

    this->postProcess(u1, this->outputNodeAttributeName1);
    this->postProcess(u2, this->outputNodeAttributeName2);
    this->mesh->setNodeAttributeFromValues("Final Cell Density Tot", 0, u1+u2);

    cout << "featk2PopulationsReactionDiffusionSolver: Info: Post processing done" << endl;
}


#endif // FEATK2POPULATIONSREACTIONDIFFUSIONSOLVER_H
