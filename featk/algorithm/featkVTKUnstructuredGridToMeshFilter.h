/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkVTKUnstructuredGridToMeshFilter.h

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
 * @class featkVTKUnstructuredGridToMeshFilter
 *
 * @brief Filter for converting a vtkUnstructuredGrid into a featkMesh.
 *
 * featkVTKUnstructuredGridToMeshFilter is a filter that produces a
 * featkMesh from a vtkUnstructuredGrid.
 *
 * featkMeshToVTKUnstructuredGridFilter translates vtkPoints into featkNode
 * coordinates, vtkCell into featkElement node connectivity, and vtkPointData
 * and vtkCellData arrays into featkNode and featkElement attributes
 * respectively.
 *
 * @warning For now, only 3D vtkUnstructuredGrid with VTK_TETRA and/or
 * VTK_HEXAHEDRON cells are supported.
 *
 * @tparam Dimension The cartesian dimension of the filter.
 *
 */

#ifndef FEATKVTKUNSTRUCTUREDGRIDTOMESHFILTER_H
#define FEATKVTKUNSTRUCTUREDGRIDTOMESHFILTER_H

#include <featk/algorithm/featkMeshProducerBase.h>
#include <featk/core/featkGlobal.h>
#include <featk/geometry/featkMesh.h>

#include <set>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

template<unsigned int Dimension>
class featkVTKUnstructuredGridToMeshFilter : public featkMeshProducerBase<3> {

    public:

        featkVTKUnstructuredGridToMeshFilter();
        ~featkVTKUnstructuredGridToMeshFilter();

        void passAllAttributesOff();
        void passAllAttributesOn();
        void setElementAttributesToPass(std::set<std::string> attributes);
        void setInputVTKUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> mesh);
        void setNodeAttributesToPass(std::set<std::string> attributes);
        void setPassAllAttributes(bool pass);
        FEATK_EXPORT void execute();

    private:

        std::set<std::string> elementAttributesToPass;
        vtkSmartPointer<vtkUnstructuredGrid> input;
        std::set<std::string> nodeAttributesToPass;
        bool passAllAttributes;
};

template<unsigned int Dimension>
featkVTKUnstructuredGridToMeshFilter<Dimension>::featkVTKUnstructuredGridToMeshFilter() {

    this->input = nullptr;
    this->passAllAttributes = true;
}

template<unsigned int Dimension>
featkVTKUnstructuredGridToMeshFilter<Dimension>::~featkVTKUnstructuredGridToMeshFilter() {

}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::passAllAttributesOff() {

    this->passAllAttributes = false;
}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::passAllAttributesOn() {

    this->passAllAttributes = true;
}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::setElementAttributesToPass(std::set<std::string> attributes) {

    this->elementAttributesToPass = attributes;
}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::setInputVTKUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> mesh) {

    this->input = mesh;
}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::setNodeAttributesToPass(std::set<std::string> attributes) {

    this->nodeAttributesToPass = attributes;
}

template<unsigned int Dimension>
void featkVTKUnstructuredGridToMeshFilter<Dimension>::setPassAllAttributes(bool pass) {

    this->passAllAttributes = pass;
}

#endif // FEATKVTKUNSTRUCTUREDGRIDTOMESHFILTER_H
