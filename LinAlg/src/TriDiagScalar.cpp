#include <iostream>
#include <iomanip>

#include "TriDiagScalar.h"  

TriDiagScalar::TriDiagScalar()  : TriDiagArray() {} 

// Konstruktor mit drei Konstanten für die Diagonalelemente
TriDiagScalar::TriDiagScalar(int lda, double diag, double upperDiag, double lowerDiag) 
	: TriDiagArray(lda), diagValue(diag), upperValue(upperDiag), lowerValue(lowerDiag) {}


// Die Matrix multipliziert mit einem Vektor
TNT::Vector<double> TriDiagScalar::ATimesX(TNT::Vector<double> x) const
{
	// Erstes Element gleich im Konstruktor besetzen
	TNT::Vector<double> result(this->n, diagValue * x[0] + upperValue * x[1]);

	for (unsigned int i = 1; i < this->n - 1; i++)
		result[i] = lowerValue * x[i - 1] + diagValue * x[i] + upperValue * x[i + 1];

	// Letztes Element wieder ohne Schleife
	result[n - 1] = lowerValue * x[n - 2] + diagValue * x[n - 1];
	return(result);
}

// Ausgabe auf stdout
void TriDiagScalar::print() const 
{
	std::cout << "Eine tridiagonale Matrix mit konstanten Werten" << std::endl;
	// Defaultanzahl Nachkommastellen abfragen, um sie
	// am Ende der Funktion wieder zu rekonstruieren
	std::streamsize defaultPrecision = std::cout.precision();
	std::cout.precision(3);

	std::cout << "Die Anzahl der Zeilen und Spalten ist ";
	std::cout << std::fixed << "n = " << this->n << std::endl;
	std::cout << "Der konstante Wert auf der Hauptdiagonale ist " << std::setprecision(3)  << diagValue << std::endl;
	std::cout << "Der konstante Wert auf der oberen Nebendiagonale ist " << std::setprecision(3) << upperValue << std::endl;
	std::cout << "Der konstante Wert auf der unteren Nebenddiagonale ist " << std::setprecision(3) << lowerValue << std::endl;

	std::cout.precision(defaultPrecision);
}


// Abfragen der oberen Nebendiagonale
TNT::Vector<double> TriDiagScalar::getDiagonal() const
{
	TNT::Vector<double> result(this->n, diagValue);
	return result;
}

// Abfragen der Elemente der unteren Nebendiagonale der Matrix
TNT::Vector<double> TriDiagScalar::getUpperDiagonal() const
{
	TNT::Vector<double> result(this->n, upperValue);
	return result;
}

// Abfragen des i-ten Elements der unteren Nebendiagonale
TNT::Vector<double> TriDiagScalar::getLowerDiagonal() const
{
	TNT::Vector<double> result(this->n, lowerValue);
	return result;
}