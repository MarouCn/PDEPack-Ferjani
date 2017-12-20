/**!
 * \file kuehlrippeTest.cpp
 * \brief Programme de test.
 * \author Marouane Ferjani
 * \version 0.1
 * \date 6 Dezember 20017
 * \headerfile hallo
 *
 * Diese Klasse ist die Loesung Programm , hier werden die 3 Algortihmus durchgefuert und dazu die VTK Datein auch fuer jede Loesung.
 *
 */

#include <cmath>
#include <iostream>
#include <C:\tnt\tnt.h>
#include "Resolve.h"

/*!
* \example kuelrippeTest.cpp
*
*
*  Diese Klasse ist die Loesung Programm , hier werden die 3 Algortihmus durchgefuert und dazu die VTK Datein auch fuer jede Loesung.
*/


int main(void)
{

        std::cout << "*************************HPC Pruefung : Loesung Programm fuer das Kuehlrippe Model***************  " << std::endl;
	
        std::cout << " " << std::endl;
        std::cout << "//////// Loesung mit dem Gauss Seidel Verfahren.//////////" << std::endl;
        Resolve* oResolve = new Resolve(90,	150, 0, false);
        oResolve->calculateGS();
        std::cout << "Loesung mit dem Gauss Seidel Verfahren beendet mit der  Dauer: " << std::fixed << oResolve->getDuration() << std::endl;

        if(!oResolve->appendToLog("LoesungBenchmark.txt")) std::cout << "Fehler beim Schreiben der Benchmark-Datei!" << std::endl;

        if (oResolve->exportToVTK("Loesung-GS.vtk"))
	{
		std::cout << "VTK-Datei erstellt!" << std::endl;
	}
	else std::cout << "Fehler beim Schreiben der VTK-Datei!" << std::endl;
        delete oResolve;

        

        std::cout << " " << std::endl;

        std::cout << "//////////Loesung mit dem Jacobi Verfahren////////////." << std::endl;
        oResolve = new Resolve(90, 150, 0, false);
        oResolve->calculateJ();
        std::cout << "Loesung mit dem Jacobi Verfahren beendet mit der Dauer: " << std::fixed << oResolve->getDuration() << std::endl;

        if (!oResolve->appendToLog("LoesungBenchmark.txt")) std::cout << "Fehler beim Schreiben der Benchmark-Datei!" << std::endl;

        if (oResolve->exportToVTK("Loesung-J.vtk"))
	{
		std::cout << "VTK-Datei erstellt!" << std::endl;
	}
	else std::cout << "Fehler beim Schreiben der VTK-Datei!" << std::endl;
        delete oResolve;

        std::cout << "Loesung mit dem Jacobi Verfahren mit ----OpenMp." << std::endl;
        oResolve = new Resolve(90, 150, 0, true);
        oResolve->calculateJ();
        std::cout << "Loesung mit dem Jacobi Verfahren mit ----OpenMp beendet mit der Dauer: " << std::fixed << oResolve->getDuration() << std::endl;

        if (!oResolve->appendToLog("LoesungBenchmark.txt")) std::cout << "Fehler beim Schreiben der Benchmark-Datei!" << std::endl;

        if (oResolve->exportToVTK("Loesung-J-openmp.vtk"))
	{
		std::cout << "VTK-Datei erstellt!" << std::endl;
	}
	else std::cout << "Fehler beim Schreiben der VTK-Datei!" << std::endl;
        delete oResolve;

        std::cout << " " << std::endl;

        std::cout << "//////Loesung mit dem SOR Verfahren ////////." << std::endl;
        oResolve = new Resolve(90, 150, 0, false);
        oResolve->calculateSOR();
        std::cout << "Loesung mit dem SOR Verfahren beendet mit der Dauer: " << std::fixed << oResolve->getDuration() << std::endl;

        if (!oResolve->appendToLog("LoesungBenchmark.txt")) std::cout << "Fehler beim Schreiben der Benchmark-Datei!" << std::endl;

        if (oResolve->exportToVTK("Loesung-sor.vtk"))
	{
		std::cout << "VTK-Datei erstellt!" << std::endl;
	}
	else std::cout << "Fehler beim Schreiben der VTK-Datei!" << std::endl;
        delete oResolve;

        std::cout << "Loesung mit dem SOR Verfahren mit ----OpenMp." << std::endl;
        oResolve = new Resolve(90, 150, 0, true);
        oResolve->calculateSOR();
        std::cout << "Loesung mit dem Jacobi Verfahren mit ----OpenMp beendet mit der Dauer: " << std::fixed << oResolve->getDuration() << std::endl;

        if (!oResolve->appendToLog("LoesungBenchmark.txt")) std::cout << "Fehler beim Schreiben der Benchmark-Datei!" << std::endl;

        if (oResolve->exportToVTK("Loesung-sor-openmp.vtk"))
	{
		std::cout << "VTK-Datei erstellt!" << std::endl;
	}
	else std::cout << "Fehler beim Schreiben der VTK-Datei!" << std::endl;
        delete oResolve;

        std::cout << " " << std::endl;

        std::cout << "End of Program!" << std::endl;

	system("pause");
	return EXIT_SUCCESS;
}
