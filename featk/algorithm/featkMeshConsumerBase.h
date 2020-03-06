/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkMeshConsumerBase.h

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
 * @class featkMeshConsumerBase
 *
 * @brief Base class for mesh consuming algorithms.
 *
 * @tparam Dimension The cartesian dimension of the algorithm.
 *
 */

#ifndef FEATKMESHCONSUMERBASE_H
#define FEATKMESHCONSUMERBASE_H

#include <featk/algorithm/featkAlgorithmBase.h>
#include <featk/geometry/featkMesh.h>

#include <vector>

template<unsigned int Dimension>
class featkMeshConsumerBase : public featkAlgorithmBase {

    public:

        virtual ~featkMeshConsumerBase();

        void setInputMesh(featkMesh<Dimension>* mesh, unsigned int port=0);

    protected:

        featkMeshConsumerBase();

        std::vector<featkMesh<Dimension>*> inputMeshes;
};

#endif // FEATKMESHCONSUMERBASE_H
