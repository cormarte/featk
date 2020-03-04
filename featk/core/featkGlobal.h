/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkGlobals.h

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
 * @brief featkGlobal defines export/import declarations for the Finite Element
 * Analysis Toolkit.
 *
 */

#ifndef FEATKGLOBAL_H
#define FEATKGLOBAL_H

#if defined(FEATK_LIBRARY)
#  define FEATK_EXPORT __declspec(dllexport)
#else
#  define FEATK_EXPORT __declspec(dllimport)
#endif

#endif // FEATKGLOBAL_H