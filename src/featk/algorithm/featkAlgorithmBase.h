/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkAlgorithmBase.h

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
 * @class featkAlgorithmBase
 *
 * @brief Base class for algorithms.
 *
 */

#ifndef FEATKALGORITHMBASE_H
#define FEATKALGORITHMBASE_H

#include <featk/core/featkGlobal.h>

class FEATK_EXPORT featkAlgorithmBase {

    public:

        virtual ~featkAlgorithmBase();

        virtual bool check();
        virtual void execute()=0;

        bool update();

    protected:

        featkAlgorithmBase();
};

#endif // FEATKALGORITHMBASE_H
