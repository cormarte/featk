#ifndef FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H
#define FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H

#include <featk/geometry/featkMesh.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

template<unsigned int Dimension>
class featkMeshToVTKUnstructuredGridFilter {

    public:

        featkMeshToVTKUnstructuredGridFilter();
        ~featkMeshToVTKUnstructuredGridFilter();

        vtkSmartPointer<vtkUnstructuredGrid> getOutputVTKUnstructuredGrid();
        void setInputMesh(featkMesh<Dimension>* mesh);
        void update();

    private:

        featkMesh<Dimension>* input;
        vtkSmartPointer<vtkUnstructuredGrid> output;
};

template<unsigned int Dimension>
featkMeshToVTKUnstructuredGridFilter<Dimension>::featkMeshToVTKUnstructuredGridFilter() {

    this->input = nullptr;
    this->output = nullptr;
}

template<unsigned int Dimension>
featkMeshToVTKUnstructuredGridFilter<Dimension>::~featkMeshToVTKUnstructuredGridFilter() {

}

template<unsigned int Dimension>
vtkSmartPointer<vtkUnstructuredGrid> featkMeshToVTKUnstructuredGridFilter<Dimension>::getOutputVTKUnstructuredGrid() {

    return this->output;
}

template<unsigned int Dimension>
void featkMeshToVTKUnstructuredGridFilter<Dimension>::setInputMesh(featkMesh<Dimension>* mesh) {

    this->input = mesh;
}

#endif // FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H
