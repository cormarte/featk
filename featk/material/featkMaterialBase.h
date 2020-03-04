#ifndef FEATKMATERIALBASE_H
#define FEATKMATERIALBASE_H

#include <featk/core/featkDefines.h>

template<unsigned int Dimension>
class featkMaterialBase {

    public:

        virtual ~featkMaterialBase();

        AttributeValueType<Dimension, 4> getConstitutiveMatrix();

    protected:

        featkMaterialBase();

        AttributeValueType<Dimension, 4> constitutiveMatrix;
};

template<unsigned int Dimension>
featkMaterialBase<Dimension>::featkMaterialBase() {

}

template<unsigned int Dimension>
featkMaterialBase<Dimension>::~featkMaterialBase() {

}

template<unsigned int Dimension>
 AttributeValueType<Dimension, 4> featkMaterialBase<Dimension>::getConstitutiveMatrix() {

    return this->constitutiveMatrix;
}

#endif // FEATKMATERIALBASE_H
