#include <featk/algorithm/featk3DGridSource.h>
#include <featk/geometry/featkHex8Element.h>
#include <featk/geometry/featkNode.h>

#include <vector>

featk3DGridSource::featk3DGridSource() {

    this->dimensions = {10, 151, 10};
    this->elementType = FEATK_HEX8;
    this->origin = {0.0, 0.0, 0.0};
    this->spacing = {1.0, 1.0, 1.0};
}

featk3DGridSource::~featk3DGridSource() {

}

void featk3DGridSource::execute() {

    // Nodes

    AttributeValueType<3, 1> origin = (AttributeValueType<3, 1>() << this->origin[0], this->origin[1], this->origin[2]).finished();

    unsigned int nx = this->dimensions[0];
    unsigned int ny = this->dimensions[1];
    unsigned int nz = this->dimensions[2];

    std::vector<featkNode<3>*> nodes;
    size_t id = 0;

    for (size_t z=0; z!=nz; z++) {

        for (size_t y=0; y!=ny; y++) {

            for (size_t x=0; x!=nx; x++) {

                AttributeValueType<3, 1> coordinates = (AttributeValueType<3, 1>() << this->spacing[0]*x, this->spacing[1]*y, this->spacing[2]*z).finished()-origin;

                featkNode<3>* node = new featkNode<3>(id, coordinates);
                nodes.push_back(node);

                id++;
            }
        }
    }


    // Elements

    std::vector<featkElementInterface<3>* > elements;

    for (size_t z=0; z!=nz-1; z++) {

        for (size_t y=0; y!=ny-1; y++) {

            for (size_t x=0; x!=nx-1; x++) {

                size_t n0 =     z*nx*ny +     y*nx +     x;
                size_t n1 =     z*nx*ny +     y*nx + (x+1);
                size_t n2 =     z*nx*ny + (y+1)*nx + (x+1);
                size_t n3 =     z*nx*ny + (y+1)*nx +     x;
                size_t n4 = (z+1)*nx*ny +     y*nx +     x;
                size_t n5 = (z+1)*nx*ny +     y*nx + (x+1);
                size_t n6 = (z+1)*nx*ny + (y+1)*nx + (x+1);
                size_t n7 = (z+1)*nx*ny + (y+1)*nx +     x;

                switch (this->elementType) {

                    case FEATK_HEX8: {

                            featkHex8Element* element = new featkHex8Element(std::vector<featkNode<3>*>({nodes[n0], nodes[n1], nodes[n2], nodes[n3], nodes[n4],  nodes[n5],  nodes[n6], nodes[n7]}));
                            elements.push_back(element);
                        }
                        break;

                case FEATK_TET4: {

                            std::vector<std::vector<featkNode<3>*>> tetrahedronNodes = {{nodes[n0], nodes[n2], nodes[n5], nodes[n1]},
                                                                                        {nodes[n2], nodes[n7], nodes[n5], nodes[n6]},
                                                                                        {nodes[n0], nodes[n7], nodes[n2], nodes[n3]},
                                                                                        {nodes[n0], nodes[n5], nodes[n2], nodes[n7]},
                                                                                        {nodes[n0], nodes[n5], nodes[n7], nodes[n4]}};

                            for (std::vector<featkNode<3>*> elementNodes : tetrahedronNodes) {

                                featkTet4Element* element = new featkTet4Element(elementNodes);
                                elements.push_back(element);
                            }
                        }
                        break;
                }
            }
        }
    }


    // Mesh

    featkMesh<3>* mesh = new featkMesh<3>(nodes, elements);


    // Output

    this->outputMeshes[0] = mesh;
}

void featk3DGridSource::setDimensions(std::array<unsigned int, 3> dimensions) {

    this->dimensions = dimensions;
}

void featk3DGridSource::setElementType(featkElementType type) {

    this->elementType = type;
}

void featk3DGridSource::setElementTypeToFEATKHex8() {

    this->elementType = FEATK_HEX8;
}

void featk3DGridSource::setElementTypeToFEATKTet4() {

    this->elementType = FEATK_TET4;
}

void featk3DGridSource::setOrigin(std::array<double, 3> origin) {

    this->origin = origin;
}

void featk3DGridSource::setSpacing(std::array<double, 3> spacing) {

    this->spacing = spacing;
}
