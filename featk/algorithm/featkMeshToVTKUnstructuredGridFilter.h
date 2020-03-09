#ifndef FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H
#define FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H

#include <featk/algorithm/featkMeshConsumerBase.h>
#include <featk/core/featkGlobal.h>
#include <featk/geometry/featkMesh.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

template<unsigned int Dimension>
class featkMeshToVTKUnstructuredGridFilter : public featkMeshConsumerBase<3> {

    public:

        featkMeshToVTKUnstructuredGridFilter();
        ~featkMeshToVTKUnstructuredGridFilter();

        vtkSmartPointer<vtkUnstructuredGrid> getOutputVTKUnstructuredGrid();
        FEATK_EXPORT void execute();

    private:

        vtkSmartPointer<vtkUnstructuredGrid> output;
};

template<unsigned int Dimension>
featkMeshToVTKUnstructuredGridFilter<Dimension>::featkMeshToVTKUnstructuredGridFilter() {

    this->output = nullptr;
}

template<unsigned int Dimension>
featkMeshToVTKUnstructuredGridFilter<Dimension>::~featkMeshToVTKUnstructuredGridFilter() {

}

template<unsigned int Dimension>
vtkSmartPointer<vtkUnstructuredGrid> featkMeshToVTKUnstructuredGridFilter<Dimension>::getOutputVTKUnstructuredGrid() {

    return this->output;
}

#endif // FEATKMESHTOVTKUNSTRUCTUREDGRIDFILTER_H
