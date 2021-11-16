#include <featk/core/featkDefines.h>
#include <featk/geometry/featkHex8Element.h>
#include <featk/geometry/featkMesh.h>
#include <featk/geometry/featkNode.h>
#include <featk/geometry/featkTet4Element.h>
#include <featk/material/featkIsotropicLinearElastic3DMaterial.h>
#include <featk/solve/featkBoundaryConditions.h>
#include <featk/solve/featkLinearElasticitySolver.h>
#include <featk/test/featkTests.h>

#include <vector>

using namespace std;

bool featkHex8StiffnessMatrixTest() {

    /**
     * From C. Felippa. Advanced Finite Element Methods: Chapter 11 - The 8-Node Hexahedron, p.16, Figure 11.28. 2017.
     */

    const unsigned int Dimension = 3;
    const unsigned int Nodes = 8;
    const unsigned int Order = 1;

    KMatrixType<Dimension, Nodes, Order> assertion = (KMatrixType<Dimension, Nodes, Order>() <<  16.0,  6.0,  6.0, -8.0,  2.0,  2.0, -6.0, -6.0,  1.0,  4.0, -2.0,  3.0,  4.0,  3.0, -2.0, -6.0,  1.0, -6.0, -4.0, -3.0, -3.0,  0.0, -1.0, -1.0,
                                                                                                  6.0, 16.0,  6.0, -2.0,  4.0,  3.0, -6.0, -6.0,  1.0,  2.0, -8.0,  2.0,  3.0,  4.0, -2.0, -1.0,  0.0, -1.0, -3.0, -4.0, -3.0,  1.0, -6.0, -6.0,
                                                                                                  6.0,  6.0, 16.0, -2.0,  3.0,  4.0, -1.0, -1.0,  0.0,  3.0, -2.0,  4.0,  2.0,  2.0, -8.0, -6.0,  1.0, -6.0, -3.0, -3.0, -4.0,  1.0, -6.0, -6.0,
                                                                                                 -8.0, -2.0, -2.0, 16.0, -6.0, -6.0,  4.0,  2.0, -3.0, -6.0,  6.0, -1.0, -6.0, -1.0,  6.0,  4.0, -3.0,  2.0,  0.0,  1.0,  1.0, -4.0,  3.0,  3.0,
                                                                                                  2.0,  4.0,  3.0, -6.0, 16.0,  6.0, -2.0, -8.0,  2.0,  6.0, -6.0,  1.0,  1.0,  0.0, -1.0, -3.0,  4.0, -2.0, -1.0, -6.0, -6.0,  3.0, -4.0, -3.0,
                                                                                                  2.0,  3.0,  4.0, -6.0,  6.0, 16.0, -3.0, -2.0,  4.0,  1.0, -1.0,  0.0,  6.0,  1.0, -6.0, -2.0,  2.0, -8.0, -1.0, -6.0, -6.0,  3.0, -3.0, -4.0,
                                                                                                 -6.0, -6.0, -1.0,  4.0, -2.0, -3.0, 16.0,  6.0, -6.0, -8.0,  2.0, -2.0, -4.0, -3.0,  3.0,  0.0, -1.0,  1.0,  4.0,  3.0,  2.0, -6.0,  1.0,  6.0,
                                                                                                 -6.0, -6.0, -1.0,  2.0, -8.0, -2.0,  6.0, 16.0, -6.0, -2.0,  4.0, -3.0, -3.0, -4.0,  3.0,  1.0, -6.0,  6.0,  3.0,  4.0,  2.0, -1.0,  0.0,  1.0,
                                                                                                  1.0,  1.0,  0.0, -3.0,  2.0,  4.0, -6.0, -6.0, 16.0,  2.0, -3.0,  4.0,  3.0,  3.0, -4.0, -1.0,  6.0, -6.0, -2.0, -2.0, -8.0,  6.0, -1.0, -6.0,
                                                                                                  4.0,  2.0,  3.0, -6.0,  6.0,  1.0, -8.0, -2.0,  2.0, 16.0, -6.0,  6.0,  0.0,  1.0, -1.0, -4.0,  3.0, -3.0, -6.0, -1.0, -6.0,  4.0, -3.0, -2.0,
                                                                                                 -2.0, -8.0, -2.0,  6.0, -6.0, -1.0,  2.0,  4.0, -3.0, -6.0, 16.0, -6.0, -1.0, -6.0,  6.0,  3.0, -4.0,  3.0,  1.0,  0.0,  1.0, -3.0,  4.0,  2.0,
                                                                                                  3.0,  2.0,  4.0, -1.0,  1.0,  0.0, -2.0, -3.0,  4.0,  6.0, -6.0, 16.0,  1.0,  6.0, -6.0, -3.0,  3.0, -4.0, -6.0, -1.0, -6.0,  2.0, -2.0, -8.0,
                                                                                                  4.0,  3.0,  2.0, -6.0,  1.0,  6.0, -4.0, -3.0,  3.0,  0.0, -1.0,  1.0, 16.0,  6.0, -6.0, -8.0,  2.0, -2.0, -6.0, -6.0, -1.0,  4.0, -2.0, -3.0,
                                                                                                  3.0,  4.0,  2.0, -1.0,  0.0,  1.0, -3.0, -4.0,  3.0,  1.0, -6.0,  6.0,  6.0, 16.0, -6.0, -2.0,  4.0, -3.0, -6.0, -6.0, -1.0,  2.0, -8.0, -2.0,
                                                                                                 -2.0, -2.0, -8.0,  6.0, -1.0, -6.0,  3.0,  3.0, -4.0, -1.0,  6.0, -6.0, -6.0, -6.0, 16.0,  2.0, -3.0,  4.0,  1.0,  1.0,  0.0, -3.0,  2.0,  4.0,
                                                                                                 -6.0, -1.0, -6.0,  4.0, -3.0, -2.0,  0.0,  1.0, -1.0, -4.0,  3.0, -3.0, -8.0, -2.0,  2.0, 16.0, -6.0,  6.0,  4.0,  2.0,  3.0, -6.0,  6.0,  1.0,
                                                                                                  1.0,  0.0,  1.0, -3.0,  4.0,  2.0, -1.0, -6.0,  6.0,  3.0, -4.0,  3.0,  2.0,  4.0, -3.0, -6.0, 16.0, -6.0, -2.0, -8.0, -2.0,  6.0, -6.0, -1.0,
                                                                                                 -6.0, -1.0, -6.0,  2.0, -2.0, -8.0,  1.0,  6.0, -6.0, -3.0,  3.0, -4.0, -2.0, -3.0,  4.0,  6.0, -6.0, 16.0,  3.0,  2.0,  4.0, -1.0,  1.0,  0.0,
                                                                                                 -4.0, -3.0, -3.0,  0.0, -1.0, -1.0,  4.0,  3.0, -2.0, -6.0,  1.0, -6.0, -6.0, -6.0,  1.0,  4.0, -2.0,  3.0, 16.0,  6.0,  6.0, -8.0,  2.0,  2.0,
                                                                                                 -3.0, -4.0, -3.0,  1.0, -6.0, -6.0,  3.0,  4.0, -2.0, -1.0,  0.0, -1.0, -6.0, -6.0,  1.0,  2.0, -8.0,  2.0,  6.0, 16.0,  6.0, -2.0,  4.0,  3.0,
                                                                                                 -3.0, -3.0, -4.0,  1.0, -6.0, -6.0,  2.0,  2.0, -8.0, -6.0,  1.0, -6.0, -1.0, -1.0,  0.0,  3.0, -2.0,  4.0,  6.0,  6.0, 16.0, -2.0,  3.0,  4.0,
                                                                                                  0.0,  1.0,  1.0, -4.0,  3.0,  3.0, -6.0, -1.0,  6.0,  4.0, -3.0,  2.0,  4.0,  2.0, -3.0, -6.0,  6.0, -1.0, -8.0, -2.0, -2.0, 16.0, -6.0, -6.0,
                                                                                                 -1.0, -6.0, -6.0,  3.0, -4.0, -3.0,  1.0,  0.0, -1.0, -3.0,  4.0, -2.0, -2.0, -8.0,  2.0,  6.0, -6.0,  1.0,  2.0,  4.0,  3.0, -6.0, 16.0,  6.0,
                                                                                                 -1.0, -6.0, -6.0,  3.0, -3.0, -4.0,  6.0,  1.0, -6.0, -2.0,  2.0, -8.0, -3.0, -2.0,  4.0,  1.0, -1.0,  0.0,  2.0,  3.0,  4.0, -6.0,  6.0, 16.0).finished();

    vector<featkNode<3>*> nodes = {new featkNode<3>(0, (AttributeValueType<Dimension, 1>() << -1.0, -1.0, -1.0).finished()),
                                   new featkNode<3>(1, (AttributeValueType<Dimension, 1>() <<  1.0, -1.0, -1.0).finished()),
                                   new featkNode<3>(2, (AttributeValueType<Dimension, 1>() <<  1.0,  1.0, -1.0).finished()),
                                   new featkNode<3>(3, (AttributeValueType<Dimension, 1>() << -1.0,  1.0, -1.0).finished()),
                                   new featkNode<3>(4, (AttributeValueType<Dimension, 1>() << -1.0, -1.0,  1.0).finished()),
                                   new featkNode<3>(5, (AttributeValueType<Dimension, 1>() <<  1.0, -1.0,  1.0).finished()),
                                   new featkNode<3>(6, (AttributeValueType<Dimension, 1>() <<  1.0,  1.0,  1.0).finished()),
                                   new featkNode<3>(7, (AttributeValueType<Dimension, 1>() << -1.0,  1.0,  1.0).finished())};

    vector<featkElementInterface<3>*> elements = {new featkHex8Element(nodes)};

    featkIsotropicLinearElastic3DMaterial material = featkIsotropicLinearElastic3DMaterial(32.0, 1.0/3.0);

    featkMesh<3>* mesh = new featkMesh<3>(nodes, elements);
    size_t id = mesh->setElementAttributeFromValues("Stiffness Tensor", 4, material.getConstitutiveMatrix());

    KMatrixType<Dimension, Nodes, Order> k = elements[0]->getBtCBIntegralMatrix<Order>(id);

    bool result = k.isApprox(assertion, EPS);

    return result;
}

void featkRunAllTests() {

    cout << featkHex8StiffnessMatrixTest() << endl;
    cout << featkTet4StiffnessMatrixTest() << endl;
    cout << featkTet4LinearElasticitySolverTest() << endl;
}

bool featkTet4StiffnessMatrixTest() {

    /**
     * From C. Felippa. Advanced Finite Element Methods: Chapter 09 - The Linear Tetrahedron, p.17, Figure 9.10. 2017.
     */

    const unsigned int Dimension = 3;
    const unsigned int Nodes = 4;
    const unsigned int Order = 1;

    KMatrixType<Dimension, Nodes, Order> assertion = (KMatrixType<Dimension, Nodes, Order>() <<  745.0,   540.0,  120.0,   -5.0,   30.0,   60.0, -270.0,  -240.0,    0.0, -470.0,  -330.0, -180.0,
                                                                                                 540.0,  1720.0,  270.0, -120.0,  520.0,  210.0, -120.0, -1080.0,  -60.0, -300.0, -1160.0, -420.0,
                                                                                                 120.0,   270.0,  565.0,    0.0,  150.0,  175.0,    0.0,  -120.0, -270.0, -120.0,  -300.0, -470.0,
                                                                                                  -5.0,  -120.0,    0.0,  145.0,  -90.0,  -60.0,  -90.0,   120.0,    0.0,  -50.0,    90.0,   60.0,
                                                                                                  30.0,   520.0,  150.0,  -90.0,  220.0,   90.0,   60.0,  -360.0,  -60.0,    0.0,  -380.0, -180.0,
                                                                                                  60.0,   210.0,  175.0,  -60.0,   90.0,  145.0,    0.0,  -120.0,  -90.0,    0.0,  -180.0, -230.0,
                                                                                                -270.0,  -120.0,    0.0,  -90.0,   60.0,    0.0,  180.0,     0.0,    0.0,  180.0,    60.0,    0.0,
                                                                                                -240.0, -1080.0, -120.0,  120.0, -360.0, -120.0,    0.0,   720.0,    0.0,  120.0,   720.0,  240.0,
                                                                                                   0.0,   -60.0, -270.0,    0.0,  -60.0,  -90.0,    0.0,     0.0,  180.0,    0.0,   120.0,  180.0,
                                                                                                -470.0,  -300.0, -120.0,  -50.0,    0.0,    0.0,  180.0,   120.0,    0.0,  340.0,   180.0,  120.0,
                                                                                                -330.0, -1160.0, -300.0,   90.0, -380.0, -180.0,   60.0,   720.0,  120.0,  180.0,   820.0,  360.0,
                                                                                                -180.0,  -420.0, -470.0,   60.0, -180.0, -230.0,    0.0,   240.0,  180.0,  120.0,   360.0,  520.0).finished();

    vector<featkNode<3>*> nodes = {new featkNode<3>(0, (AttributeValueType<Dimension, 1>() << 2.0, 3.0, 4.0).finished()),
                                   new featkNode<3>(1, (AttributeValueType<Dimension, 1>() << 6.0, 3.0, 2.0).finished()),
                                   new featkNode<3>(2, (AttributeValueType<Dimension, 1>() << 2.0, 5.0, 1.0).finished()),
                                   new featkNode<3>(3, (AttributeValueType<Dimension, 1>() << 4.0, 3.0, 6.0).finished())};


    vector<featkElementInterface<3>*> elements = {new featkTet4Element(nodes)};

    featkIsotropicLinearElastic3DMaterial material = featkIsotropicLinearElastic3DMaterial(480.0, 1.0/3.0);

    featkMesh<3>* mesh = new featkMesh<3>(nodes, elements);
    size_t id = mesh->setElementAttributeFromValues("Stiffness Tensor", 4, material.getConstitutiveMatrix());

    KMatrixType<Dimension, Nodes, Order> k = elements[0]->getBtCBIntegralMatrix<Order>(id);

    bool result = k.isApprox(assertion, EPS);

    return result;
}

bool featkTet4LinearElasticitySolverTest() {

    /**
     * From I. M. Smith, D. V. Griffiths and L. Margets. Programming the Finite Element Method, 5th Ed.: Chapter 05 - Static Equilibrium of Linear Elastic Solids, p.202, Figure 5.30. 2014.
     */

    const unsigned int Dimension = 3;
    const unsigned int Order = 1;

    vector<AttributeValueType<Dimension, Order>> assertion = {(AttributeValueType<Dimension, Order>() << 0.0000E+00, 0.0000E+00, -0.1000E-01).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.3000E-02, 0.0000E+00, -0.1000E-01).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.0000E+00, 0.0000E+00,  0.0000E+00).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.3000E-02, 0.0000E+00,  0.0000E+00).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.0000E+00, 0.3000E-02, -0.9999E-02).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.3000E-02, 0.3000E-02, -0.1000E-01).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.0000E+00, 0.3000E-02,  0.0000E+00).finished(),
                                                              (AttributeValueType<Dimension, Order>() << 0.3000E-02, 0.3000E-02,  0.0000E+00).finished()};

    vector<featkNode<3>*> nodes = {new featkNode<3>(0, (AttributeValueType<Dimension, 1>() << 0.0, 0.0,  0.0).finished()),
                                   new featkNode<3>(1, (AttributeValueType<Dimension, 1>() << 1.0, 0.0,  0.0).finished()),
                                   new featkNode<3>(2, (AttributeValueType<Dimension, 1>() << 0.0, 0.0, -1.0).finished()),
                                   new featkNode<3>(3, (AttributeValueType<Dimension, 1>() << 1.0, 0.0, -1.0).finished()),
                                   new featkNode<3>(4, (AttributeValueType<Dimension, 1>() << 0.0, 1.0,  0.0).finished()),
                                   new featkNode<3>(5, (AttributeValueType<Dimension, 1>() << 1.0, 1.0,  0.0).finished()),
                                   new featkNode<3>(6, (AttributeValueType<Dimension, 1>() << 0.0, 1.0, -1.0).finished()),
                                   new featkNode<3>(7, (AttributeValueType<Dimension, 1>() << 1.0, 1.0, -1.0).finished())};

    vector<vector<featkNode<3>*>> connectivity = {{nodes[0], nodes[3], nodes[2], nodes[6]},  // Remark: featkTet4Element uses counterclockwise convention for nodes
                                                  {nodes[0], nodes[1], nodes[3], nodes[6]},  // numbering whereas clockwise convention is used in the reference book,
                                                  {nodes[0], nodes[4], nodes[1], nodes[6]},  // leading to a sign inversion of the computed displacements. The second
                                                  {nodes[5], nodes[7], nodes[3], nodes[6]},  // and third nodes of each tetrahedron were thus inverted in this example
                                                  {nodes[5], nodes[3], nodes[1], nodes[6]},  // w.r.t. the book.
                                                  {nodes[5], nodes[1], nodes[4], nodes[6]}};

    vector<featkElementInterface<3>*> elements;

    for (std::vector<featkNode<3>*> elementNodes : connectivity) {

        featkTet4Element* element = new featkTet4Element(elementNodes);
        elements.push_back(element);
    }

    featkIsotropicLinearElastic3DMaterial material = featkIsotropicLinearElastic3DMaterial(100.0, 0.3);

    featkMesh<3>* mesh = new featkMesh<3>(nodes, elements);
    mesh->setElementAttributeFromValues("Stiffness Tensor", 4, material.getConstitutiveMatrix());
    mesh->setNodeAttributeFromValues("Body Force", 1, AttributeValueType<Dimension, 1>::Zero());

    featkBoundaryConditions<Dimension, Order>* essentialBoundaryConditions = new featkBoundaryConditions<Dimension, Order>();
    essentialBoundaryConditions->setDOFValue(0, 0, 0.0);
    essentialBoundaryConditions->setDOFValue(0, 1, 0.0);
    essentialBoundaryConditions->setDOFValue(1, 1, 0.0);
    essentialBoundaryConditions->setDOFValue(2, 0, 0.0);
    essentialBoundaryConditions->setDOFValue(2, 1, 0.0);
    essentialBoundaryConditions->setDOFValue(2, 2, 0.0);
    essentialBoundaryConditions->setDOFValue(3, 1, 0.0);
    essentialBoundaryConditions->setDOFValue(3, 2, 0.0);
    essentialBoundaryConditions->setDOFValue(4, 0, 0.0);
    essentialBoundaryConditions->setDOFValue(6, 0, 0.0);
    essentialBoundaryConditions->setDOFValue(6, 2, 0.0);
    essentialBoundaryConditions->setDOFValue(7, 2, 0.0);

    featkBoundaryConditions<Dimension, Order>* naturalBoundaryConditions = new featkBoundaryConditions<Dimension, Order>();

    AttributeValueType<Dimension, 1> f1 = (AttributeValueType<Dimension, 1>() << 0.0,  0.0, -0.1667).finished();
    AttributeValueType<Dimension, 1> f2 = (AttributeValueType<Dimension, 1>() << 0.0,  0.0, -0.3333).finished();

    for (unsigned int i=0; i!=3; i++) {

        naturalBoundaryConditions->setDOFValue(0, i, f1(i, 0));
        naturalBoundaryConditions->setDOFValue(1, i, f2(i, 0));
        naturalBoundaryConditions->setDOFValue(4, i, f2(i, 0));
        naturalBoundaryConditions->setDOFValue(5, i, f1(i, 0));
    }

    featkLinearElasticitySolver<Dimension> solver = featkLinearElasticitySolver<Dimension>();
    solver.setInputMesh(mesh);
    solver.setEssentialBoundaryConditions(essentialBoundaryConditions);
    solver.setNaturalBoundaryConditions(naturalBoundaryConditions);
    solver.update();

    size_t id = mesh->getNodeAttributeID("Displacements", Order);

    bool result = true;

    for (size_t i=0; i!=nodes.size(); i++) {

        AttributeValueType<Dimension, Order> displacements = nodes[i]->getAttributeValue(id);

        if (!displacements.isApprox(assertion[i], EPS)) {

            result = false;
        }
    }

    return result;
}
