#include <iomanip>

#include "utils.h"
#include "JacobiTriScalar.h"  

JacobiTriScalar::JacobiTriScalar() : Iterative() {}

// Konstruktor mit Koeffizientenmatrix
JacobiTriScalar::JacobiTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs) 
	: Iterative(rhs), diag(matrix.getDiagonalElement()), upper(matrix.getUpperElement()),
	  lower(matrix.getLowerElement()), xalt(this->n) 
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	if (fabs(diag) < Eps)
			std::cerr << "Das Diagonallement im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

JacobiTriScalar::JacobiTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs, double epsilon, unsigned int maxI)
	: Iterative(rhs, epsilon, maxI), diag(matrix.getDiagonalElement()), upper(matrix.getUpperElement()),
	lower(matrix.getLowerElement()), xalt(this->n)
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	if (fabs(diag) < Eps)
		std::cerr << "Das Diagonallement im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

// Wir überprüfen nicht, ob die Dimension von Start zum linearen Gleichungssystem passt
void JacobiTriScalar::iterate()
{
	xalt = x;

	// Erstes Element separat
	x[0] = (b[0] - upper*xalt[1]) / diag;
	// Jetzt die Elemente von 1 bis n-2
	for (unsigned int i = 1; i < n - 1; i++)
		x[i] = (b[i] - lower*xalt[i - 1] - upper*xalt[i + 1]) / diag;
	// Letztes Element separat
	x[n - 1] = (b[n-1] - lower*xalt[n - 2]) / diag;
	
	if (compareVectorMaximumNorm(xalt, x, localEpsilon))
		stop();
}

void JacobiTriScalar::print() const
{
	std::cout << "Gesamtschrittverfahren fuer lineare Gleichungssysteme" << std::endl;
	std::cout << "Die Koeffizientenmatrix ist eine tridiagonale Matrix mit konstanten Elementen." << std::endl;

	std::cout << "Die Koeffizientenmatrix" << std::endl;
	std::cout << "Das Diagonalelement " << diag << std::endl;
	std::cout << "Das Element auf der oberen Nebendiagonale " << upper << std::endl;
	std::cout << "Das Element auf der unteren Nebendiagonale " << lower << std::endl;

	std::cout << std::endl << "Die rechte Seite" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << b;

	Iterative::print();
}