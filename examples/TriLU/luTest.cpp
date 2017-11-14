// ---------------------------------------------------------------------------
// Einlesen und Lösen eines tridiagonalen LGS
// ---------------------------------------------------------------------------
#include <cmath>
#include <iostream>

#include <tnt.h>

#include "TriDiagArray.h"
#include "TriLU.h"

using namespace std;

/*!
 * \example luTest.cpp
 *
 * Beispiel 
 */
int main(void)
{
	int n(5);

	TriDiagArray tri(n, 2.0, -1.0, -1.0);

	// LU instanziieren und Lösung ausgeben
	TriLU lu(tri);

	// LU-Zerlegung ausgeben
	lu.print();

	// Rechte Seite
	TNT::Vector<double> rhs(n, 0.0);
	rhs[0] = 1.0;
	TNT::Vector<double> x = lu.solve(rhs);

	std::cout << "Die berechnete Loesung" << std::endl;
	std::cout << x;

	// Probe
	TNT::Vector<double> correct(tri.ATimesX(x));
	std::cout << "Probe!" << std::endl;
	std::cout << tri.ATimesX(x);

	TNT::Vector<double> a(n, 4.0), b(n - 1, -1.0), c(n - 1, -1.0);
	TriDiagArray tri2(a, b, c);

	//tri2.print();

    return 0;
}
