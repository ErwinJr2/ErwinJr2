#ifndef __SCIENCE
#define __SCIENCE
#include <stdint.h>
#define sq(X) ((X)*(X)) /**< \f$ X^2 \f$ */

#define hbar 1.0545718e-34 /**< \f$ 1.0545718 \times 10^{-34} J \cdot s \f$ */
#define m0 9.10938356e-31  /**< \f$ 9.10938356 \times 10^{-31} \f$ kg */
#define e0 1.60217662e-19  /**< \f$ 1.60217662 \times 10^{-19} \f$ C */
#define eps0 8.854187817620e-22 /**< \f$ 8.854187817620 \times 10^{-22} \f$ F/Angstrom */
#define kb 8.6173303e-5 /**< \f$ 8.6173303 \times 10^{-5} \f$ eV/K */
/* #define pi 3.1415926535897932385 --> M_PI in math.h*/
#define ANG 1E-10  /**< \f$ 10^{-10} \f$ m */
/* Angstrom in meter: all unit for length is Angstrom in the program */

typedef int32_t numpyint;
#ifdef _WINDLL
#define M_PI 3.14159265358979323846
#endif
#endif /* ifndef __SCIENCE */
