#ifndef FEATKVTKUNSTRUCTUREDGRIDTOMESHFILTER_H
#define FEATKVTKUNSTRUCTUREDGRIDTOMESHFILTER_H

#include <featk/geometry/featkMesh.h>

#include <set>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

template<unsigned int Dimension>
class featkVTKUnstructuredGridToMeshFilter {

    public:

        featkVTKUnstructuredGridToMeshFilter();
        ~featkVTKUnstructuredGridToMeshFilter();

        featkMesh<Dimension>* getOutputMesh();
        void passAllAttributesOff();
        void passAllAttributesOn();
        void setElementAttributesToPass(std::set<std::string> attributes);
        void setInputVTKUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> mesh);
        void setNodeAttributesToPass(std::set<std::string> attributes);
        void setPassAllAttributes(bool pass);
        void update();

    private:

        std::set<std::string> elementAttributesToPass;
        vtkSmartPointer<vtkUnstructuredGrid> input;
        std::set<std::string> nodeAttributesToPass;
        featkMesh<Dimension>* output;
        bool passAllAttributes;
};

template<unsigned int Dimension>
featkVTKUnstructuredGridToMeshFilter<Dimension>::featkVTKUnstructuredGridToMeshFilter() {

    this->input = nullptr;
    this->output = nullptr;
    this->passAllAttributes = true;
}

template<unsigned int Dimension>
featkVTKUnstructuredGridToMeshFilter<Dimension>::~featkVTKUnstructuredGridToMeshFilter() {

}

template<unsigned int Dimension>
featkMesh<Dimension>* featkVTKUnstructuredGridToMeshFilter<Dimension>::getOutputMesh() {

    return this->output;
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
