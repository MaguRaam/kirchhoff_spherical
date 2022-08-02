#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include "Headers.h"

void ThrowOrderError(unsigned int); 
void ThrowNegativePressureDensityError(Vector<double>, Vector<double>, Point<2>, bool); 


#endif /* EXCEPTIONS_H_ */
