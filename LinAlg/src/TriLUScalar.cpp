#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "utils.h"
#include "TriLUScalar.h"  

TriLUScalar::TriLUScalar() : TriLU() {}

TriLUScalar::TriLUScalar(unsigned int lda, double diag, double upper, double lower)
	: TriLU(lda), M(lda, diag, upper, lower) {}


// Konstruktor mit Koeffizientenmatrix
TriLUScalar::TriLUScalar(TriDiagScalar matrix)
	: TriLU(matrix.getDimension()), M(matrix)
{
	this->u = TNT::Vector<double>(n-1, matrix.getUpperElement());
	short e = decompose();
	if (e != 0)
		std::cerr << "Fehler bei der Berechnung der LU-Zerlegung in TriLUScalar: Diagonalelement zu klein!" << std::endl;
}

// LU-Zerlegung mit konstanten Werten in den Diagonalen
short TriLUScalar::decompose()
{
	// Die Werte auf den Diagonalen abfragen
	double a(M.getDiagonalElement()),
		   b(M.getUpperElement()),
		   c(M.getLowerElement());


	d[0] = a;

	for (unsigned int i = 0; i < this->n - 1; i++) {
		if (fabs(d[i]) > Eps)
			l[i] = c / d[i];
		else {
			return(1);
		}
		d[i + 1] = a - l[i] * b;
	}
	return(0);
}

void TriLUScalar::print() const
{
	std::cout << "LU-Zerlegung eines tridiagonalen linearen Gleichungssystems " << std::endl;
	// Defaultanzahl Nachkommastellen abfragen, um sie
	// am Ende der Funktion wieder zu rekonstruieren
	std::streamsize defaultPrecision = std::cout.precision();
	std::cout.precision(3);

	std::cout << "Die Koeffizientenmatrix" << std::endl;
	M.print();

	std::cout << std::endl << "Die untere Diagonale von L" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << l;
	std::cout << "Die Diagonale von U" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << d;

	std::cout << std::endl <<  "Die obere Diagonale von U" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << u;

	std::cout.precision(defaultPrecision);
}