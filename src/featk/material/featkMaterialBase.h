/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkMaterialBase.h

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
 * @class featkMaterialBase
 *
 * @brief Base class for materials.
 *
 * featkMaterialBase is a base class for materials. featkMaterialBase stores
 * the material 4th-order constitutive tensor used in continuous mechanics
 * finite elements problems.
 *
 * @tparam Dimension The cartesian dimension of the material.
 *
 */

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
