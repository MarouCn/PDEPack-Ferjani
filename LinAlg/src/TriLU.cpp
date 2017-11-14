#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "utils.h"
#include "TriLU.h"  

TriLU::TriLU() : n(0) {}

TriLU::TriLU(unsigned int lda) : n(lda), d(n, 0.0), l(n - 1, 0.0), u(n - 1, 0.0) {}

// Konstruktor mit Koeffizientenmatrix
TriLU::TriLU(TriDiagArray matrix) 
	: n(matrix.getDimension()), A(matrix), d(n, 0.0), l(n - 1, 0.0), u(n-1, matrix.getUpperDiagonal())
{
	short e = decompose();
	if (e != 0)
		std::cerr << "Fehler bei der Berechnung der LU-Zerlegung in TriLU: Diagonalelement zu klein!" << std::endl;
}

// LU-Zerlegung
short TriLU::decompose()
{
	TNT::Vector<double> a(A.getDiagonal()),
		                b(A.getUpperDiagonal()),
		                c(A.getLowerDiagonal());

	d[0] = a[0];

	for (unsigned int i = 0; i < this->n - 1; i++) {
		if (fabs(d[i]) > Eps)
			l[i] = c[i] / d[i];
		else {
			return(1);
		}
		d[i + 1] = a[i + 1] - l[i] * b[i];
	}
	return(0);
}

// Die Vorwärtssubstitution
TNT::Vector<double> TriLU::forwardSubstitution(TNT::Vector<double> rhs) const
{
	TNT::Vector<double> y(n, rhs[0]);

	for (unsigned int i = 1; i < this->n; i++)
		y[i] = rhs[i] - l[i - 1] * y[i - 1];

	return(y);
}

// Die Rückwärtssubstitution
TNT::Vector<double> TriLU::backwardSubstition(TNT::Vector<double> rhs) const
{
	TNT::Vector<double> x(n, 0.0);
	x[n - 1] = rhs[n-1]/ d[n - 1];
	for (int i = n - 2; i >= 0; i--)
		x[i] = (rhs[i] - u[i] * x[i + 1])/d[i];
	return x;
}

TNT::Vector<double> TriLU::solve(TNT::Vector<double> rhs) const
{
	return this->backwardSubstition(this->forwardSubstitution(rhs));
}

void TriLU::print() const
{
	std::cout << "LU-Zerlegung eines tridiagonalen linearen Gleichungssystems" << std::endl;
	// Defaultanzahl Nachkommastellen abfragen, um sie
	// am Ende der Funktion wieder zu rekonstruieren
	std::streamsize defaultPrecision = std::cout.precision();
	std::cout.precision(3);

	std::cout << "Die Koeffizientenmatrix" << std::endl;
	A.print();

	std::cout << std::endl << "Die untere Diagonale von L" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << l;
	std::cout << "Die Diagonale von U" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << d;

	std::cout << std::endl <<  "Die obere Diagonale von U" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << u;

	std::cout.precision(defaultPrecision);
}