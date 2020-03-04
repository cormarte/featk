/**
  * See Flaherty. Finite Element Analysis, ch.4, p.24. (http://www.cs.rpi.edu/~flaherje/pdf/fea4.pdf) for paramertization.
  */

#ifndef FEATKTET4ELEMENT_H
#define FEATKTET4ELEMENT_H

#include <featk/geometry/featkElement.h>

using featkTet4Element = featkElement<3, 4, 4, 3>;

#endif // FEATKTET4ELEMENT_H
