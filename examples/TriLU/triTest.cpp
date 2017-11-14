// ---------------------------------------------------------------------------
// Einlesen und Lösen eines tridiagonalen LGS
// ---------------------------------------------------------------------------
#include <cmath>
#include <iostream>

#include <tnt.h>

#include "TriDiagArray.h"
#include "TriDiagScalar.h"

using namespace std;

/*!
 * \example triTest.cpp
 *
 * Matrizen der beiden Klassen besetzen, ausgeben und die Funktion ATimesX
 * ausfuehren. Das exakte Ergebnis wird auf der Konsole ausgegeben, um
 * die korrekte Ausführung zu überprüfen.
 */
int main(void)
{
	int n(5);

	TriDiagArray tri(n, 2.0, -1.0, -1.0);

	tri.print();

	// ATimesX testen
	TNT::Vector<double> x(n, 1.0);
	std::cout << "Test der Multiplikation mit einem Vektor" << std::endl;
	std::cout << "Das Ergebnis muss der Vektor (1,0,0,0,1) sein!" << std::endl;
	std::cout << tri.ATimesX(x);

	TNT::Vector<double> a(n, 4.0), b(n - 1, -1.0), c(n - 1, -1.0);
	TriDiagArray tri2(a, b, c);

	tri2.print();

	// ATimesX testen
	std::cout << "Test der Multiplikation mit einem Vektor" << std::endl;
	std::cout << "Das Ergebnis muss der Vektor (3, 2, 2, 2,  3) sein!" << std::endl;
	std::cout << tri2.ATimesX(x);


	// Und jetzt noch die Variante, die nur Konstante in der Instanz speichert
	TriDiagScalar tri3(n, 4.0, -1.0, 1.0);
	tri3.print();
	// ATimesX testen
	std::cout << "Test der Multiplikation mit einem Vektor" << std::endl;
	std::cout << "Das Ergebnis muss der Vektor (3, 4, 4, 4,  5) sein!" << std::endl;
	std::cout << tri3.ATimesX(x);

    return 0;
}
