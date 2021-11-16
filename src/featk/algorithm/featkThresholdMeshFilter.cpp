#include <featk/algorithm/featkThresholdMeshFilter.h>

#include <iostream>

void featkThresholdMeshFilter<3>::execute() {

    featkMesh<3>* inputMesh = this->inputMeshes[0];

    if (inputMesh != nullptr) {

        size_t attributeID = inputMesh->getElementAttributeID(this->attributeName, 0);

        if (attributeID != 0) {

            // Thresholding

            std::vector<featkElementInterface<3>*> lowerElements;
            std::vector<featkElementInterface<3>*> upperElements;
            std::vector<featkNode<3>*> lowerNodes;
            std::vector<featkNode<3>*> upperNodes;
            this->nodeMaps = std::vector<std::map<size_t, size_t>>(2);
            this->elementMaps = std::vector<std::map<size_t, size_t>>(2);

            for (size_t e=0; e!=inputMesh->getNumberOfElements(); e++) {

                featkElementInterface<3>* inputElement = inputMesh->getElement(e);

                double value = inputElement->getAttributeValue(attributeID)(0,0);

                std::vector<featkElementInterface<3>*>& elements = value < this->threshold ? lowerElements : upperElements;
                std::vector<featkNode<3>*>& nodes = value < this->threshold ? lowerNodes : upperNodes;
                std::map<size_t, size_t>& elementMap = value < this->threshold ? this->elementMaps[0] : this->elementMaps[1];
                std::map<size_t, size_t>& nodeMap = value < this->threshold ? this->nodeMaps[0] : this->nodeMaps[1];

                std::vector<featkNode<3>*> elementNodes;

                for (featkNode<3>* inputNode : inputElement->getNodes()) {

                    size_t id = inputNode->getID();
                    featkNode<3>* node;

                    if (nodeMap.count(id)){

                        node = nodes[nodeMap[id]];
                    }

                    else {

                        node = new featkNode<3>(nodes.size(), inputNode->getCoordinates());
                        nodeMap[id] = node->getID();
                        nodes.push_back(node);
                    }

                    elementNodes.push_back(node);
                }

                featkElementInterface<3>* element;

                switch (inputElement->getElementType()) {

                    case FEATK_TET4:

                        element = new featkTet4Element(elementNodes);
                        break;

                    case FEATK_HEX8:

                        element = new featkHex8Element(elementNodes);
                        break;

                    default:

                        element = nullptr;
                        std::cout << "featkThresholdMeshFilter: Error: Element type not yet supported." << std::endl;
                        break;
                }

                elementMap[e] = elements.size();
                elements.push_back(element);
            }

            this->outputMeshes[0] = new featkMesh<3>(lowerNodes, lowerElements);
            this->outputMeshes[1] = new featkMesh<3>(upperNodes, upperElements);


            // Reversed maps

            this->reversedNodeMaps = std::vector<std::map<size_t, size_t>>(2);

            for (int i=0; i!=this->reversedNodeMaps.size(); i++) {

                for (auto it=this->nodeMaps[i].begin(); it!=this->nodeMaps[i].end(); ++it) {

                    this->reversedNodeMaps[i][it->second] = it->first;
                }
            }


            this->reversedElementMaps = std::vector<std::map<size_t, size_t>>(2);

            for (int i=0; i!=this->reversedElementMaps.size(); i++) {

                for (auto it=this->elementMaps[i].begin(); it!=this->elementMaps[i].end(); ++it) {

                    this->reversedElementMaps[i][it->second] = it->first;
                }
            }


            // Copying attributes

            for (const auto& pair : inputMesh->getNodeAttributeTable()) {

                std::string name = pair.first;
                size_t id = pair.second.first;
                unsigned int order = pair.second.second;

                for (int i=0; i!=this->outputMeshes.size(); i++) {

                    featkMesh<3>* outputMesh = this->outputMeshes[i];

                    std::vector<std::shared_ptr<MatrixXd>> attributes;

                    for (size_t j=0; j!=outputMesh->getNumberOfNodes(); j++) {

                        attributes.push_back(std::make_shared<MatrixXd>(inputMesh->getNode(this->reversedNodeMaps[i][j])->getAttributeValue(id)));
                    }

                    outputMesh->setNodeAttributes(name, order, attributes);
                }
            }

            for (const auto& pair : inputMesh->getElementAttributeTable()) {

                std::string name = pair.first;
                size_t id = pair.second.first;
                unsigned int order = pair.second.second;

                for (int i=0; i!=this->outputMeshes.size(); i++) {

                    featkMesh<3>* outputMesh = this->outputMeshes[i];

                    std::vector<std::shared_ptr<MatrixXd>> attributes;

                    for (size_t j=0; j!=outputMesh->getNumberOfElements(); j++) {

                        attributes.push_back(std::make_shared<MatrixXd>(inputMesh->getElement(this->reversedElementMaps[i][j])->getAttributeValue(id)));
                    }

                    outputMesh->setElementAttributes(name, order, attributes);
                }
            }
        }
    }
}
