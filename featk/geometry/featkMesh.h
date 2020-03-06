#ifndef FEATKMESH_H
#define FEATKMESH_H

#include <featk/geometry/featkElementInterface.h>
#include <featk/geometry/featkNode.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

template<unsigned int Dimension>
class featkMesh {

    public:

        featkMesh(std::vector<featkNode<Dimension>*> nodes, std::vector<featkElementInterface<Dimension>*> elements);
        ~featkMesh();

        template<unsigned int Order> void computeNodeBQ(std::string inputNodeAttributeName, std::string outputNodeAttributeName);
        template<unsigned int Order> void computeNodeCBQ(std::string elementAttributeName, std::string inputNodeAttributeName, std::string outputNodeAttributeName);

        void addNodeAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values);
        featkElementInterface<Dimension>* getElement(size_t index) const;
        size_t getElementAttributeID(std::string name, unsigned int order) const;
        std::map<std::string, std::pair<size_t, unsigned int>> getElementAttributeTable() const;
        MatrixXd getElementAttributeValues(std::string name, unsigned int order) const;
        const std::vector<featkElementInterface<Dimension>*>& getElements() const;
        featkNode<Dimension>* getNode(size_t index) const;
        size_t getNodeAttributeID(std::string name, unsigned int order) const;
        std::map<std::string, std::pair<size_t, unsigned int>> getNodeAttributeTable() const;
        MatrixXd getNodeAttributeValues(std::string name, unsigned int order) const;
        const std::vector<featkNode<Dimension>*>& getNodes() const;
        size_t getNumberOfElements() const;
        size_t getNumberOfNodes() const;
        void setElementAttributes(std::string name, unsigned int order, const std::vector<std::shared_ptr<MatrixXd>>& attributes);
        size_t setElementAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values);
        void setNodeAttributes(std::string name, unsigned int order, const std::vector<std::shared_ptr<MatrixXd>>& attributes);
        size_t setNodeAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values);

    private:

        template<typename AttributeOwnerType> void addAttributes(const std::vector<AttributeOwnerType*>& items, size_t id, std::vector<std::shared_ptr<MatrixXd>> attributes) const;
        template<typename AttributeOwnerType> MatrixXd getAttributeValues(const std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, const std::vector<AttributeOwnerType*>& items, std::string name, unsigned int order) const;
        template<typename AttributeOwnerType> void setAttributes(const std::vector<AttributeOwnerType*>& items, size_t id, std::vector<std::shared_ptr<MatrixXd>> attributes) const;

        size_t getAttributeID(const std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, std::string name, unsigned int order) const;
        std::vector<std::shared_ptr<MatrixXd>> getAttributesFromValues(size_t items, unsigned int order, const MatrixXd& values) const;
        size_t registerAttribute(std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, size_t& attributeMaxID, std::string name, unsigned int order);

        size_t elementAttributeMaxID;
        std::map<std::string, std::pair<size_t, unsigned int>> elementAttributeTable;
        std::vector<featkElementInterface<Dimension>*> elements;
        size_t nodeAttributeMaxID;
        std::map<std::string, std::pair<size_t, unsigned int>> nodeAttributeTable;
        std::vector<featkNode<Dimension>*> nodes;
};

template<unsigned int Dimension>
featkMesh<Dimension>::featkMesh(std::vector<featkNode<Dimension>*> nodes, std::vector<featkElementInterface<Dimension>*> elements) {

    this->elements = elements;
    this->nodes = nodes;

    this->elementAttributeMaxID = 0;
    this->nodeAttributeMaxID = 0;

    this->registerAttribute(this->nodeAttributeTable, this->nodeAttributeMaxID, "Cartesian Coordinates", 1);
}

template<unsigned int Dimension>
featkMesh<Dimension>::~featkMesh() {

    for (auto* element : this->elements) {

        delete element;
    }

    for (auto* node : this->nodes) {

        delete node;
    }
}


template<unsigned int Dimension>
template<typename AttributeOwnerType>
void featkMesh<Dimension>::addAttributes(const std::vector<AttributeOwnerType*>& items, size_t id, std::vector<std::shared_ptr<MatrixXd>> attributes) const {

    if (items.size() == attributes.size()) {

        for (size_t i=0; i!=items.size(); i++) {

            items[i]->addAttribute(id, attributes[i]);
        }
    }
}

template<unsigned int Dimension>
template<unsigned int Order>
void featkMesh<Dimension>::computeNodeBQ(std::string inputNodeAttributeName, std::string outputNodeAttributeName) {

    MatrixXd values;

    unsigned int rows = POWER(Dimension, (Order+1));
    unsigned int cols = 1;

    size_t id = this->getNodeAttributeID(inputNodeAttributeName, Order);

    if (id != 0) {

        values = MatrixXd(this->nodes.size()*rows, cols);

        size_t index = 0;

        for (featkNode<Dimension>* node : this->nodes) {

            MatrixXd value = MatrixXd::Zero(rows, cols);
            std::vector<featkElementInterface<Dimension>*> elements = node->getElements();

            for (featkElementInterface<Dimension>* element : elements) {

                value += element->getNodeBQVector<Order>(node, id);
            }

            value /= elements.size();
            values.block(index, 0, rows, cols) = value;

            index += rows;
        }
    }

    else {

        values = MatrixXd::Zero(this->nodes.size()*rows, cols);
    }

    this->setNodeAttributeFromValues(outputNodeAttributeName, Order+1, values);
}

template<unsigned int Dimension>
template<unsigned int Order>
void featkMesh<Dimension>::computeNodeCBQ(std::string elementAttributeName, std::string inputNodeAttributeName, std::string outputNodeAttributeName) {

    MatrixXd values;

    unsigned int rows = POWER(Dimension, (Order+1));
    unsigned int cols = 1;

    size_t elementAttributeID = this->getElementAttributeID(elementAttributeName, 2*(Order+1));
    size_t inputNodeAttributeID = this->getNodeAttributeID(inputNodeAttributeName, Order);

    if (elementAttributeID != 0 && inputNodeAttributeID != 0) {

        values = MatrixXd(this->nodes.size()*rows, cols);

        size_t index = 0;

        for (featkNode<Dimension>* node : this->nodes) {

            MatrixXd value = MatrixXd::Zero(rows, cols);
            std::vector<featkElementInterface<Dimension>*> elements = node->getElements();

            for (featkElementInterface<Dimension>* element : elements) {

                value += element->getNodeCBQVector<Order>(node, elementAttributeID, inputNodeAttributeID);
            }

            value /= elements.size();
            values.block(index, 0, rows, cols) = value;

            index += rows;
        }
    }

    else {

        values = MatrixXd::Zero(this->nodes.size()*rows, cols);
    }

    this->setNodeAttributeFromValues(outputNodeAttributeName, Order+1, values);
}

template<unsigned int Dimension>
template<typename AttributeOwnerType>
MatrixXd featkMesh<Dimension>::getAttributeValues(const std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, const std::vector<AttributeOwnerType*>& items, std::string name, unsigned int order) const {

    MatrixXd values;

    unsigned int rows = POWER(Dimension, order/2+order%2);
    unsigned int cols = POWER(Dimension, order/2);

    size_t id = this->getAttributeID(attributeTable, name, order);

    if (id != 0) {

        values = MatrixXd(items.size()*rows, cols);

        size_t index = 0;

        for (const auto& item : items) {

            values.block(index, 0, rows, cols) = item->getAttributeValue(id);
            index += rows;
        }
    }

    else {

        values = MatrixXd::Zero(items.size()*rows, cols);
    }

    return values;
}

template<unsigned int Dimension>
template<typename AttributeOwnerType>
void featkMesh<Dimension>::setAttributes(const std::vector<AttributeOwnerType*>& items, size_t id, std::vector<std::shared_ptr<MatrixXd>> attributes) const {

    if (items.size() == attributes.size()) {

        for (size_t i=0; i!=items.size(); i++) {

            items[i]->setAttribute(id, attributes[i]);
        }
    }
}


template<unsigned int Dimension>
void featkMesh<Dimension>::addNodeAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values) {

    size_t id = this->getNodeAttributeID(name, order);

    if (id != 0) {

        std::vector<std::shared_ptr<MatrixXd>> attributes = this->getAttributesFromValues(this->nodes.size(), order, values);  // This must be checked before registering
        this->addAttributes(this->nodes, id, attributes);
    }
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::getAttributeID(const std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, std::string name, unsigned int order) const {

    return attributeTable.count(name) && attributeTable.at(name).second == order ? attributeTable.at(name).first : 0;
}

template<unsigned int Dimension>
std::vector<std::shared_ptr<MatrixXd>> featkMesh<Dimension>::getAttributesFromValues(size_t items, unsigned int order, const MatrixXd& values) const {

    unsigned int components = POWER(Dimension, order);
    unsigned int rows = POWER(Dimension, order/2+order%2);
    unsigned int cols = POWER(Dimension, order/2);

    std::vector<std::shared_ptr<MatrixXd>> attributes(items);

    if (values.cols()==1) {

        if (values.rows()==components*items) {

            size_t index = 0;

            for (size_t i=0; i!=items; i++) {

                MatrixXd value = MatrixXd(rows, cols);

                for (unsigned int j=0; j!=rows; j++) {

                    for (unsigned int k=0; k!=cols; k++) {

                        value(j, k) = values(index, 0);
                        index++;
                    }
                }

                attributes[i] = std::make_shared<MatrixXd>(value);
            }
        }

        else if (values.rows()==components) {

            MatrixXd value = MatrixXd(rows, cols);

            size_t index = 0;

            for (unsigned int j=0; j!=rows; j++) {

                for (unsigned int k=0; k!=cols; k++) {

                    value(j, k) = values(index, 0);
                    index++;
                }
            }

            for (size_t i=0; i!=items; i++) {

                attributes[i] = std::make_shared<MatrixXd>(value);
            }
        }

        else {

            /**
             * Degenerated case
             */
        }
    }

    else if (values.cols()==cols) {

        if (values.rows()==rows*items) {

            size_t index = 0;

            for (size_t i=0; i!=items; i++) {

                MatrixXd value = MatrixXd(values.block(index, 0, rows, cols));
                index+= rows;

                attributes[i] = std::make_shared<MatrixXd>(value);
            }
        }

        else if (values.rows()==rows) {

            for (size_t i=0; i!=items; i++) {

                attributes[i] = std::make_shared<MatrixXd>(values);
            }
        }

        else {

            /**
             * Degenerated case
             */
        }
    }

    else {

        /**
         * Degenerated case
         */
    }

    return attributes;
}

template<unsigned int Dimension>
featkElementInterface<Dimension>* featkMesh<Dimension>::getElement(size_t index) const {

    return this->elements[index];
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::getElementAttributeID(std::string name, unsigned int order) const {

    return this->getAttributeID(this->elementAttributeTable, name, order);
}

template<unsigned int Dimension>
std::map<std::string, std::pair<size_t, unsigned int>> featkMesh<Dimension>::getElementAttributeTable() const {

    return this->elementAttributeTable;
}

template<unsigned int Dimension>
MatrixXd featkMesh<Dimension>::getElementAttributeValues(std::string name, unsigned int order) const {

    return this->getAttributeValues(this->elementAttributeTable, this->elements, name, order);
}

template<unsigned int Dimension>
const std::vector<featkElementInterface<Dimension>*>& featkMesh<Dimension>::getElements() const {

    return this->elements;
}

template<unsigned int Dimension>
featkNode<Dimension>* featkMesh<Dimension>::getNode(size_t index) const {

    return this->nodes[index];
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::getNodeAttributeID(std::string name, unsigned int order) const {

    return this->getAttributeID(this->nodeAttributeTable, name, order);
}

template<unsigned int Dimension>
std::map<std::string, std::pair<size_t, unsigned int>> featkMesh<Dimension>::getNodeAttributeTable() const {

    return this->nodeAttributeTable;
}

template<unsigned int Dimension>
MatrixXd featkMesh<Dimension>::getNodeAttributeValues(std::string name, unsigned int order) const {

    return this->getAttributeValues(this->nodeAttributeTable, this->nodes, name, order);
}

template<unsigned int Dimension>
const std::vector<featkNode<Dimension>*>& featkMesh<Dimension>::getNodes() const {

    return this->nodes;
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::getNumberOfElements() const {

    return this->elements.size();
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::getNumberOfNodes() const {

    return this->nodes.size();
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::registerAttribute(std::map<std::string, std::pair<size_t, unsigned int>>& attributeTable, size_t& attributeMaxID, std::string name, unsigned int order) {

    if (attributeTable.count(name)) {

        attributeTable[name].second = order;
    }

    else {

        attributeTable[name] = std::make_pair(++attributeMaxID, order);
    }

    return attributeTable[name].first;
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::setElementAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values) {

    std::vector<std::shared_ptr<MatrixXd>> attributes = this->getAttributesFromValues(this->elements.size(), order, values);  // This must be checked before registering
    size_t id = this->registerAttribute(this->elementAttributeTable, this->elementAttributeMaxID, name, order);
    this->setAttributes(this->elements, id, attributes);

    return id;
}

template<unsigned int Dimension>
size_t featkMesh<Dimension>::setNodeAttributeFromValues(std::string name, unsigned int order, const MatrixXd& values) {

    std::vector<std::shared_ptr<MatrixXd>> attributes = this->getAttributesFromValues(this->nodes.size(), order, values);  // This must be checked before registering
    size_t id = this->registerAttribute(this->nodeAttributeTable, this->nodeAttributeMaxID, name, order);
    this->setAttributes(this->nodes, id, attributes);

    return id;
}

#endif // FEATKMESH_H
