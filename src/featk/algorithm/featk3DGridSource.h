/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featk3DGridSource.h

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
 * @class featk3DGridSource
 *
 * @brief 3D grid mesh source.
 *
 */

#ifndef FEATK3DGRIDSOURCE_H
#define FEATK3DGRIDSOURCE_H

#include <featk/algorithm/featkMeshProducerBase.h>
#include <featk/core/featkGlobal.h>
#include <featk/geometry/featkMesh.h>

#include <array>

class FEATK_EXPORT featk3DGridSource : public featkMeshProducerBase<3> {

    public:

        featk3DGridSource();
        ~featk3DGridSource();

        void execute();

        void setDimensions(std::array<unsigned int, 3> dimensions);
        void setElementType(featkElementType type);
        void setElementTypeToFEATKHex8();
        void setElementTypeToFEATKTet4();
        void setOrigin(std::array<double, 3> origin);
        void setSpacing(std::array<double, 3> spacing);

    private:

        std::array<unsigned int, 3> dimensions;
        featkElementType elementType;
        std::array<double, 3> origin;
        std::array<double, 3> spacing;
};

#endif // FEATK3DGRIDSOURCE_H
