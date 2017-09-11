// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef GOLDEN_SECTION_SEARCH_H
#define GOLDEN_SECTION_SEARCH_H
 
// This is the content of the .h file, which is where the declarations go
#include <valarray>
using std::valarray;

valarray<double> golden_section_search(int n, valarray<double> x, double a, valarray<double> d); 
 
// This is the end of the header guard
#endif
