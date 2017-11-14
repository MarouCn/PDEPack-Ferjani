#include <iomanip>

#include "utils.h"
#include "Jacobi.h"  

Jacobi::Jacobi() : Iterative() {}

// Konstruktor mit Koeffizientenmatrix
Jacobi::Jacobi(TNT::Matrix<double> matrix, TNT::Vector<double> rhs) 
	: Iterative(rhs), A(matrix), xalt(n) 
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	for (unsigned int i = 0; i < this->n; i++)
		if (fabs(A[i][i]) < Eps)
			std::cerr << "Das Diagonallement " << i << " im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

Jacobi::Jacobi(TNT::Matrix<double> matrix, TNT::Vector<double> rhs, double epsilon, unsigned int maxI)
	: Iterative(rhs, epsilon, maxI), A(matrix), xalt(n) 
{
	// Wir überprüfen, ob die Diagonalelemente nicht zu klein sind
	for (unsigned int i = 0; i < this->n; i++)
		if (fabs(A[i][i]) < Eps)
			std::cerr << "Das Diagonallement " << i << " im Gesamtschritt-Verfahren ist zu klein!" << std::endl;
}

// Wir überprüfen nicht, ob die Dimension von Start zum linearen Gleichungssystem passt
void Jacobi::iterate()
{
	unsigned int i, j;
	double sum(0.0);

	// Werte von x auf xalt kopieren
	xalt = x;
	for (i = 0; i < n; i++)
	{
		sum = 0.0;
		for (j = 0; j < i; j++)
			sum += A[i][j] * xalt[j];
		for (j = i + 1; j < n; j++)
			sum += A[i][j] * xalt[j];
		// Wir haben im Konstruktor überprüft, dass die
		// Diagonalelemente ok sind.
		x[i] = (b[i] - sum) / A[i][i];

	}
	
	if (compareVectorMaximumNorm(xalt, x, localEpsilon))
		stop();
}

void Jacobi::print() const
{
	std::cout << "Gesamtschrittverfahren fuer lineare Gleichungssysteme" << std::endl;

	std::cout << "Die Koeffizientenmatrix" << std::endl;
	std::cout << A << std::endl;

	std::cout << std::endl << "Die rechte Seite" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << b;

	Iterative::print();
}