// ---------------------------------------------------------------------------
// Lösung der eindimensionalen Drahtgleichung
// Literatur: van Kan, Segal, Vermolen: 
// Numerical Methods in Scientific Computing, Absatz 3.1 "The cable equation"
// ---------------------------------------------------------------------------
#include <cmath>
#include <iostream>

#include <tnt.h>
#include "utils.h"
#include "visUtils.h"

#include "Wire1D.h"

//! \addtogroup apps Drahtgleichung
/*! \addtogroup apps
*  @{
* Lösung der Drahtgleichung für einen im Quelltext vorgegebenen Kraftvektor
* und Export der Ergebnisse im VTK-Format für die Visualisierung in ParaView.
*/

int main(void)
{
    int n(10);
	double len(11.0);
	Wire1D wire(len, n);

	// Kraftwerte
	TNT::Vector<double> f(n);
	f[0] = 0.1;
	f[1] = 0.2;
	f[2] = 0.3;
	f[3] = 0.4;
	f[4] = 0.5;
	f[5] = 0.5;
	f[6] = 0.4;
	f[7] = 0.3;
	f[8] = 0.2;
	f[9] = 0.1;

	TNT::Vector<double> y = wire.solve(f);

	// Die numersiche Lösung ausgeben, um die Paraview-Pipeline zu füllen
	storeLine("first.vtk", wire.getXCoord(), y);

    return EXIT_SUCCESS;
}
/*! @}  end doxygen addtogroup apps */