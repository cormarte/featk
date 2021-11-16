/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkMeshProducerBase.h

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
 * @class featkMeshToVTKUnstructuredGridFilter
 *
 * @brief Filter for converting a featkMesh into a vtkUnstructuredGrid.
 *
 * featkMeshToVTKUnstructuredGridFilter is a filter that produces a
 * vtkUnstructuredGrid from a featkMesh.
 *
 * featkMeshToVTKUnstructuredGridFilter translates featkNode coordinates
 * into vtkPoints, featkElement node connectivity into vtkCell, and
 * featkNode and featkElement attributes into vtkDoubleArray which are
 * respectively assigned to the vtkPointData and vtkCellData of the output
 * vtkUnstructuredGrid.
 *
 * @warning For now, only 3D featkMesh with featkTet4Element and/or
 * featkHex8Element are supported.
 *
 * @tparam Dimension The cartesian dimension of the filter.
 *
 */

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
class featkMeshToVTKUnstructuredGridFilter : public featkMeshConsumerBase<Dimension> {

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
