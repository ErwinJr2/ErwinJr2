#ifndef __SCIENCE
#define __SCIENCE 
#define sq(X) ((X)*(X))

#define hbar 1.0545718e-34 /*J.s*/
#define m0 9.10938356e-31  /*kg*/
#define e0 1.60217662e-19  /*C*/
#define eps0 8.854187817620e-22 /*F/Angstrom*/
#define kb 8.6173303e-5 /*eV/K*/
/* #define pi 3.1415926535897932385 --> M_PI in math.h*/
#define ANG 1E-10  
/* Angstrom in meter: all unit for length is Angstrom in the program */

#ifdef _WINDLL
typedef int32_t numpyint;
#else
typedef int64_t numpyint;
#endif
#endif /* ifndef __SCIENCE */
