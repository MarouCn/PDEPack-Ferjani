#ifndef RESOLVE_H
#define RESOLVE_H

#include <fstream>
#include <iostream>
#include <C:\tnt\tnt.h>
#include <string>
#include "getWalltime.h"
#include "utils.h"


//! Klasse zur Berechnung der W�rmeleitung in einer K�hlrippe. Die konstante Temperatur <temperature> wirkt von links auf die K�hlrippe ein. Die Umgebungstemperatur betr�gt 0 Grad.
class Resolve {

protected:

	//! Die Koeffizientenmatrix.
	TNT::Matrix<double> matrix;

	//! Die Dimension der Matrix
	int dimension;

	//! Die konstante Temperatur, die an die K�hlrippe abgegeben wird
	int temperature;

	//! Der konstante Funktionswert f(xi,xj)
	double f;

	//! h = 1/dimension
	double h;

	//! Boolescher Wert, ob OpenMP zur Parallelisierung verwendet werden soll
	bool useOpenMP;

	//! Die maximale Dauer der Berechnung
	int timeLimit;

	//! Die tats�chliche Dauer der Berechnung
	double duration;

	//! Die Anzahl der Schleifeniterationen
	int count;

	//! Das K�rzel des gew�hlten Verfahrens
	std::string algo = "None";

public:

	//! Der Defaultkonstruktor
        Resolve();


	
	//! Ein Konstruktor mit der einwirkenden Temperatur, der Dimension der Matrix, einem Zeitlimit f�r die Berechnung in Sekunden und einem Flag, ob OpenMP zur Parallelisierung verwendet werden soll
	/*!
	* \param temperature
	*  Die konstante Temperatur, die an die K�hlrippe abgegeben wird
	* \param dimension
	*  Die Dimension der Matrix
	* \param timeLimit
	*  Zeitlimit f�r die Berechnung in Sekunden
	* \param useOpenMP
	*  Booelscher Wert, ob OpenMP zur Parallelisierung verwendet werden soll
	*/
        Resolve(int temperature, int dimension, int timeLimit, bool useOpenMP);

	

	/*!
	* Berechnet den Temperaturverlauf unter Verwendung des SOR-Verfahrens mit dem optimalen Relaxationsparameter woptimal. Das Ergebnis wird in der Membervariable matrix gespeichert.
	*/
	void calculateSOR();

	/*!
	* Berechnet den Temperaturverlauf unter Verwendung des Gauss-Seidel-Verfahrens. Das Ergebnis wird in der Membervariable matrix gespeichert.
	*/
	void calculateGS();

	/*!
	* Berechnet den Temperaturverlauf unter Verwendung des Jacobi-Verfahrens. Das Ergebnis wird in der Membervariable matrix gespeichert.
	*/
	void calculateJ();

	/*!
	* Exportiert die Matrix in eine VTK-Datei mit dem Namen <fileName>
	*/
	bool exportToVTK(std::string fileName) const;

	/*!
	* Schreibt die Berechnungsparameter und die Berechnungsdauer in eine Log-Datei mit dem Namen <fileName>
	*/
	bool appendToLog(std::string fileName) const;

	//! Gibt die Anzahl der Schleifendurchl�ufe zur�ck.
	int getIterations();

	//! Gibt die tats�chliche Berechnungsdauer zur�ck.
	double getDuration();

	//! Gibt die Matrix zur�ck.
	TNT::Matrix<double> getMatrix();
};

#endif
