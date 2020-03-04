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
