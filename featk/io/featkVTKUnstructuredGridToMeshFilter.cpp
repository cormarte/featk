#include <featk/core/featkDefines.h>
#include <featk/geometry/featkHex8Element.h>
#include <featk/geometry/featkTet4Element.h>
#include <featk/io/featkVTKUnstructuredGridToMeshFilter.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

void featkVTKUnstructuredGridToMeshFilter<3>::update() {

    if (this->input != nullptr) {

        // Points and cells

        vtkPoints* points = this->input->GetPoints();
        vtkIdType numberOfPoints = points->GetNumberOfPoints();

        vtkCellArray* cells = this->input->GetCells();
        vtkIdType numberOfCells = cells->GetNumberOfCells();


        // Identify used points

        std::vector<bool> used(numberOfPoints);

        vtkIdType cellLocation = 0;

        for (vtkIdType i=0; i!=numberOfCells; i++) {

            vtkIdType numIds;
            vtkIdType* pointIds;
            cells->GetCell(cellLocation, numIds, pointIds);

            for (vtkIdType j=0; j!=numIds; j++) {

                used[pointIds[j]] = true;
            }

            cellLocation += 1+numIds;
        }


        // Nodes

        std::vector<featkNode<3>*> nodes;
        std::map<vtkIdType, featkNode<3>*> map;
        size_t id = 0;

        for (vtkIdType n=0; n!=numberOfPoints; n++) {

            if (used[n]) {

                double point[3];
                points->GetPoint(n, point);

                featkNode<3>* node = new featkNode<3>(id, (AttributeValueType<3, 1>() << point[0], point[1], point[2]).finished());
                nodes.push_back(node);
                map[n] = node;
                id++;
            }
        }

        size_t numberOfUsedPoints = nodes.size();
        cout << "featkVTKUnstructuredGridToMeshFilter: " << (numberOfPoints-numberOfUsedPoints) << "/" << numberOfPoints << " points unused." << endl;


        // Elements

        std::vector<featkElementInterface<3>*> elements(numberOfCells);

        cellLocation = 0;

        for (vtkIdType e=0; e!=numberOfCells; e++) {

            vtkIdType numIds;
            vtkIdType* pointIds;
            cells->GetCell(cellLocation, numIds, pointIds);

            std::vector<featkNode<3>*> elementNodes;

            for (vtkIdType j=0; j!=numIds; j++) {

                elementNodes.push_back(map[pointIds[j]]);
            }

            featkElementInterface<3>* element;

            switch (this->input->GetCellType(e)) {

                case VTK_TETRA: {

                    element = new featkTet4Element(elementNodes);
                    break;
                }

                case VTK_HEXAHEDRON: {

                    element = new featkHex8Element(elementNodes);
                    break;
                }
            }

            elements[e] = element;

            cellLocation += 1+numIds;
        }


        // Mesh

        featkMesh<3>* mesh = new featkMesh<3>(nodes, elements);


        // Node attributes

        vtkPointData* pointData = this->input->GetPointData();

        for (int i=0; i!=pointData->GetNumberOfArrays(); i++) {

            vtkDataArray* array = pointData->GetArray(i);
            std::string name = array->GetName();

            if (this->passAllAttributes || this->nodeAttributesToPass.count(name)) {

                int numberOfComponents = array->GetNumberOfComponents();

                if (IS_POWER(numberOfComponents, 3)) {

                    unsigned int order = LOG(numberOfComponents, 3);
                    VectorXd values(numberOfUsedPoints*numberOfComponents, 1);
                    id = 0;

                    for (vtkIdType j=0; j!=numberOfPoints; j++) {

                        if (used[j]) {

                            for (int k=0; k!=numberOfComponents; k++) {

                                values[id] = array->GetComponent(j, k);
                                id++;
                            }
                        }
                    }

                    mesh->setNodeAttributeFromValues(name, order, values);
                }
            }
        }


        // Element attributes

        vtkCellData* cellData = this->input->GetCellData();

        for (int i=0; i!=cellData->GetNumberOfArrays(); i++) {

            vtkDataArray* array = cellData->GetArray(i);
            std::string name = array->GetName();

            if (this->passAllAttributes || find(this->elementAttributesToPass.begin(), this->elementAttributesToPass.end(), name) != this->elementAttributesToPass.end()) {

                int numberOfComponents = array->GetNumberOfComponents();

                if (IS_POWER(numberOfComponents, 3)) {

                    unsigned int order = LOG(numberOfComponents, 3);
                    VectorXd values(numberOfCells*numberOfComponents);
                    id = 0;

                    for (vtkIdType j=0; j!=numberOfCells; j++) {

                        for (int k=0; k!=numberOfComponents; k++) {

                            values[id] = array->GetComponent(j, k);
                            id++;
                        }
                    }

                    mesh->setElementAttributeFromValues(name, order, values);
                }
            }
        }


        // Output

        this->output = mesh;
    }
}
