// ---------------------------------------------------------------------------
// Beispielprogramm für die Verwendung von Jacpbo
// mit einer tridiagonalen Koeffizientenmatrix
// ---------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>

#include <tnt.h>
#include "utils.h"

#include "JacobiTriScalar.h"

using namespace std;

/*!
* \example triJacobiTest.cpp 
*
* Einlesen eines linearen Gleichungssystems mit einer tridiagonalen Matrix
* mit konstanten Werten auf den Diagonalen und Lösung
* mit dem Jacobi-Verfahren.
*/
int main(void)
{
	ifstream input;
	string fileName;

	// Einlesen der Matrix
	fileName = "matrix3.txt";
	input.open(fileName.c_str(), std::ios::in);
	if (input.fail()) {
		std::cerr << "Eingabedatei " << fileName << " existiert nicht!" << endl;
		return 2;
	}

	TNT::Vector<double> b;
	TNT::Vector<double> start;

	input >> b;
	input >> start;
	input.close();
	TriDiagScalar A(3, 2.0, -1.0, -1.0);

	JacobiTriScalar gesamtschritt(A, b, 0.005, 1000);

	TNT::Vector<double> x = gesamtschritt.solve(start);

	gesamtschritt.print();


	std::cout << "Die berechnete Loesung" << std::endl;
	std::cout << std::fixed << std::setprecision(6) << x << std::endl;

	TNT::Vector<double> correct(b.size(), 1.0);
	std::cout << "Die exakte Loesung" << std::endl;
	std::cout << std::fixed << std::setprecision(6) << correct << std::endl;


	std::cout << std::scientific << "Die Maximumsnorm der Abweichung " << normMaximum(x - correct) << std::endl;

	return EXIT_SUCCESS;
}
