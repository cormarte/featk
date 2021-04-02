/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkMeshProducerBase.h

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
 * @class featkMeshProducerBase
 *
 * @brief Base class for mesh producing algorithms.
 *
 * featkMeshProducerBase is a base class for featkMesh producing
 * featkAlgorithm objects.
 *
 * featkMeshProducerBase stores featkMesh objects produced by the
 * featkAlgorithm.
 *
 * @tparam Dimension The cartesian dimension of the algorithm.
 *
 */

#ifndef FEATKMESHPRODUCERBASE_H
#define FEATKMESHPRODUCERBASE_H

#include <featk/algorithm/featkAlgorithmBase.h>
#include <featk/geometry/featkMesh.h>

#include <vector>

template<unsigned int Dimension>
class featkMeshProducerBase : public virtual featkAlgorithmBase {

    public:

        virtual ~featkMeshProducerBase();

        featkMesh<Dimension>* getOutputMesh(unsigned int port=0) const;

    protected:

        featkMeshProducerBase();

        void setNumberOfOutputMeshes(unsigned int number);

        std::vector<featkMesh<Dimension>*> outputMeshes;
};

template<unsigned int Dimension>
featkMeshProducerBase<Dimension>::featkMeshProducerBase() {

    this->setNumberOfOutputMeshes(1);
}

template<unsigned int Dimension>
featkMeshProducerBase<Dimension>::~featkMeshProducerBase() {

}

template<unsigned int Dimension>
void featkMeshProducerBase<Dimension>::setNumberOfOutputMeshes(unsigned int number) {

    this->outputMeshes = std::vector<featkMesh<Dimension>*>(number);
}

template<unsigned int Dimension>
featkMesh<Dimension>* featkMeshProducerBase<Dimension>::getOutputMesh(unsigned int port=0) const {

    return port < this->outputMeshes.size() ? this->outputMeshes[port] : nullptr;
}

#endif // FEATKMESHPRODUCERBASE_H
