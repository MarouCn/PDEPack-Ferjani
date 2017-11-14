#include <iostream>
#include <iomanip>

#include "TriDiagArray.h"  

TriDiagArray::TriDiagArray()  : n(0) {} 

TriDiagArray::TriDiagArray(unsigned int lda) : n(lda) {}

// Konstruktor mit drei Feldern für die Diagonalelemente
TriDiagArray::TriDiagArray(TNT::Vector<double> diag, TNT::Vector<double> upperDiag, TNT::Vector<double> lowerDiag)
	: n(diag.dim()), a(diag), b(upperDiag), c(lowerDiag) {}


// Konstruktor mit drei Konstanten für die Diagonalelemente
TriDiagArray::TriDiagArray(unsigned int lda, double diag, double upperDiag, double lowerDiag) 
	: n(lda), a(lda, diag), b(lda-1, upperDiag), c(lda-1, lowerDiag) {}


// Die Matrix multipliziert mit einem Vektor
TNT::Vector<double> TriDiagArray::ATimesX(TNT::Vector<double> x) const
{
	// Erstes Element gleich im Konstruktor besetzen
	TNT::Vector<double> result(n, a[0] * x[0] + b[0] * x[1]);
	
	for (unsigned int i = 1; i < this->n - 1; i++)
		result[i] = c[i - 1] * x[i - 1] + a[i] * x[i] + b[i] * x[i + 1];
	result[n - 1] = c[n - 2] * x[n - 2] + a[n - 1] * x[n - 1];
	return(result);
}

// Ausgabe auf stdout
void TriDiagArray::print() const 
{
	std::cout << "Eine tridiagonale Matrix" << std::endl;
	// Defaultanzahl Nachkommastellen abfragen, um sie
	// am Ende der Funktion wieder zu rekonstruieren
	std::streamsize defaultPrecision = std::cout.precision();
	std::cout.precision(3);

	std::cout << "Die Anzahl der Zeilen und Spalten ist ";
	std::cout << std::fixed << "n = " << this->n << std::endl;
	std::cout << "Die Diagonale" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << a;

	std::cout << "Die obere Diagonale" << std::endl;
	std::cout << b;

	std::cout <<  "Die untere Diagonale" << std::endl;
	std::cout << std::fixed << std::setprecision(3) << c;
	std::cout << std::endl;

	std::cout.precision(defaultPrecision);
}


// Abfragen der oberen Nebendiagonale
TNT::Vector<double> TriDiagArray::getDiagonal() const
{
	return a;
}

// Abfragen der Elemente der unteren Nebendiagonale der Matrix
TNT::Vector<double> TriDiagArray::getUpperDiagonal() const
{
	return b;
}

// Abfragen des i-ten Elements der unteren Nebendiagonale
TNT::Vector<double> TriDiagArray::getLowerDiagonal() const
{
	return c;
}