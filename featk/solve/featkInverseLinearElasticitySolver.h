#ifndef FEATKINVERSELINEARELASTICITYSOLVER_H
#define FEATKINVERSELINEARELASTICITYSOLVER_H

#include <featk/solve/featkSolverBase.h>

template<unsigned int Dimension>
class featkInverseLinearElasticitySolver : public featkSolverBase<Dimension, 1> {

    public:

        featkInverseLinearElasticitySolver();
        ~featkInverseLinearElasticitySolver();

        void solve();

    protected:

        SparseMatrix<double> getGlobalSystemMatrix();
        VectorXd getGlobalSystemVector();
        void postProcess(const VectorXd& u);

        std::string displacementAttributeName;  // Check if attributes are valid and assign their IDs to vars
        std::string outputAttributeName;
        std::string stiffnessAttributeName;
};

template<unsigned int Dimension>
featkInverseLinearElasticitySolver<Dimension>::featkInverseLinearElasticitySolver() {

    this->displacementAttributeName = "Displacements";
    this->outputAttributeName = "Force";
    this->stiffnessAttributeName = "Stiffness Tensor";
}

template<unsigned int Dimension>
featkInverseLinearElasticitySolver<Dimension>::~featkInverseLinearElasticitySolver() {

}

template<unsigned int Dimension>
SparseMatrix<double> featkInverseLinearElasticitySolver<Dimension>::getGlobalSystemMatrix() {

    return this->getGlobalMatrixFromElements(&featkInverseLinearElasticitySolver<Dimension>::getElementBtCBIntegralMatrix, {this->mesh->getElementAttributeID(this->stiffnessAttributeName, 4)});
}

template<unsigned int Dimension>
VectorXd featkInverseLinearElasticitySolver<Dimension>::getGlobalSystemVector() {

    MatrixXd q = this->mesh->getNodeAttributeValues(this->displacementAttributeName, 1);
    cout << q.rows() << ", " << q.cols() << endl;

    return q;
}

template<unsigned int Dimension>
void featkInverseLinearElasticitySolver<Dimension>::postProcess(const VectorXd& u) {

    this->mesh->setNodeAttributeFromValues(this->outputAttributeName, 1, u);
    this->mesh->computeNodeBQ<1>(this->outputAttributeName, this->outputAttributeName + " Gradient");
}

template<unsigned int Dimension>
void featkInverseLinearElasticitySolver<Dimension>::solve() {

    SparseMatrix<double> k = this->getGlobalSystemMatrix();
    VectorXd q = this->getGlobalSystemVector();
    VectorXd f = k*q;

    this->postProcess(f);
}

#endif // FEATKINVERSELINEARELASTICITYSOLVER_H
