#include <iomanip>

#include "utils.h"
#include "GaussSeidelTriScalar.h"  

GaussSeidelTriScalar::GaussSeidelTriScalar() : Iterative() {}

// Konstruktor mit Koeffizientenmatrix
GaussSeidelTriScalar::GaussSeidelTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs) 
	: Iterative(rhs), diag(matrix.getDiagonalElement()), upper(matrix.getUpperElement()),
	  lower(matrix.getLowerElement()), xalt(this->n) 
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	if (fabs(diag) < Eps)
			std::cerr << "Das Diagonallement im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

GaussSeidelTriScalar::GaussSeidelTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs, double epsilon, unsigned int maxI)
	: Iterative(rhs, epsilon, maxI), diag(matrix.getDiagonalElement()), upper(matrix.getUpperElement()),
	lower(matrix.getLowerElement()), xalt(this->n)
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	if (fabs(diag) < Eps)
		std::cerr << "Das Diagonallement im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

// Wir überprüfen nicht, ob die Dimension von Start zum linearen Gleichungssystem passt
void GaussSeidelTriScalar::iterate()
{
	xalt = x;

	// Erstes Element separat
	x[0] = (b[0] - upper*x[1]) / diag;
	// Jetzt die Elemente von 1 bis n-2
	for (unsigned int i = 1; i < n - 1; i++)
		x[i] = (b[i] - lower*x[i - 1] - upper*x[i + 1]) / diag;
	// Letztes Element separat
	x[n - 1] = (b[n-1] - lower*x[n - 2]) / diag;
	
	if (compareVectorMaximumNorm(xalt, x, localEpsilon))
		stop();
}

void GaussSeidelTriScalar::print() const
{
	std::cout << "Einzelschrittverfahren fuer lineare Gleichungssysteme" << std::endl;
	std::cout << "Die Koeffizientenmatrix ist eine tridiagonale Matrix mit konstanten Elementen." << std::endl;

	std::cout << "Die Koeffizientenmatrix" << std::endl;
	std::cout << "Das Diagonalelement " << diag << std::endl;
	std::cout << "Das Element auf der oberen Nebendiagonale " << upper << std::endl;
	std::cout << "Das Element auf der unteren Nebendiagonale " << lower << std::endl;

	std::cout << std::endl << "Die rechte Seite" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << b;

	Iterative::print();
}