#include <featk/algorithm/featkMeshToVTKUnstructuredGridFilter.h>

void featkMeshToVTKUnstructuredGridFilter<3>::execute() {

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    featkMesh<3>* inputMesh = this->inputMeshes[0];

    if (inputMesh != nullptr) {


        // Nodes allocation

        std::vector<featkNode<3>*> nodes = inputMesh->getNodes();

        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        points->SetNumberOfPoints(nodes.size());


        // Node attributes allocation

        std::vector<std::pair<size_t, vtkSmartPointer<vtkDoubleArray>>> nodeAttributes;

        for (const auto& pair : inputMesh->getNodeAttributeTable()) {

            std::string name = pair.first;
            size_t id = pair.second.first;
            unsigned int order = pair.second.second;
            unsigned int components = POWER(3, order);

            vtkSmartPointer<vtkDoubleArray> attributeArray = vtkDoubleArray::New();
            attributeArray->SetName(name.c_str());
            attributeArray->SetNumberOfComponents(components);
            attributeArray->SetNumberOfTuples(nodes.size());
            attributeArray->SetNumberOfValues(nodes.size()*components);

            nodeAttributes.push_back(std::make_pair(id, attributeArray));
        }


        // Nodes loop

        std::map<size_t, vtkIdType> map;  // featkNode<Dimension>::id is supposed to match node index in featkMesh<Dimension>::nodes but keep a map in case

        for (vtkIdType n=0; n!=nodes.size(); n++) {

            featkNode<3>* node = nodes[n];


            // Insertion

            AttributeValueType<3, 1> coordinates = node->getCoordinates();
            points->SetPoint(n, coordinates(0, 0), coordinates(1, 0), coordinates(2, 0));
            map[node->getID()] = n;


            // Attributes

            for (const auto& pair : nodeAttributes) {

                MatrixXd matrix = node->getAttributeValue(pair.first);

                double* values = new double[pair.second->GetNumberOfComponents()];
                int index = 0;

                for (int i=0; i!=matrix.rows(); i++) {

                    for (int j=0; j!=matrix.cols(); j++) {

                        values[index] = matrix(i, j);
                        index++;
                    }
                }

                pair.second->SetTuple(n, values);
                delete[] values;
            }
        }

        unstructuredGrid->SetPoints(points);


        // Node attributes assginment

        for (const auto& pair : nodeAttributes) {

            unstructuredGrid->GetPointData()->AddArray(pair.second);
        }


        // Elements allocation

        std::vector<featkElementInterface<3>*> elements = inputMesh->getElements();

        unstructuredGrid->Allocate(elements.size());


        // Element attributes allocation

        std::vector<std::pair<size_t, vtkSmartPointer<vtkDoubleArray>>> elementAttributes;

        for (const auto& pair : inputMesh->getElementAttributeTable()) {

            std::string name = pair.first;
            size_t id = pair.second.first;
            unsigned int order = pair.second.second;
            unsigned int components = POWER(3, order);

            vtkSmartPointer<vtkDoubleArray> attributeArray = vtkDoubleArray::New();
            attributeArray->SetName(name.c_str());
            attributeArray->SetNumberOfComponents(components);
            attributeArray->SetNumberOfTuples(elements.size());
            attributeArray->SetNumberOfValues(elements.size()*components);

            elementAttributes.push_back(std::make_pair(id, attributeArray));
        }


        // Elements loop

        for (vtkIdType e=0; e!=elements.size(); e++) {

            featkElementInterface<3>* element = elements[e];


            // Cell type

            int type;

            switch (element->getElementType()) {

                case FEATK_TET4:

                    type = VTK_TETRA;
                    break;

                case FEATK_HEX8:

                    type = VTK_HEXAHEDRON;
                    break;
            }


            // Insertion

            std::vector<featkNode<3>*> elementNodes = element->getNodes();

            vtkSmartPointer<vtkIdList> pointIDs = vtkIdList::New();
            pointIDs->SetNumberOfIds(elementNodes.size());

            for (vtkIdType n=0; n!=elementNodes.size(); n++) {

                pointIDs->SetId(n, map[elementNodes[n]->getID()]);
            }

            unstructuredGrid->InsertNextCell(type, pointIDs);


            // Attributes

            for (const auto& pair : elementAttributes) {

                MatrixXd matrix = element->getAttributeValue(pair.first);

                double* values = new double[pair.second->GetNumberOfComponents()];
                int index = 0;

                for (int i=0; i!=matrix.rows(); i++) {

                    for (int j=0; j!=matrix.cols(); j++) {

                        values[index] = matrix(i, j);
                        index++;
                    }
                }

                pair.second->SetTuple(e, values);
                delete[] values;
            }
        }


        // Element attributes assginment

        for (const auto& pair : elementAttributes) {

            unstructuredGrid->GetCellData()->AddArray(pair.second);
        }
    }


    // Output

    this->output = unstructuredGrid;
}
