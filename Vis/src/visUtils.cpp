#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include "visUtils.h"

short storeLine(string fileName, TNT::Vector<double> x, TNT::Vector<double> y) 
{
	int i, n(x.size());
	double z(0.0);
	ofstream output;

	// open ostream
	output.open(fileName.c_str(), std::ios::out);
	if (output.fail()) {
		std::cerr << "Fehler beim Öffnen der Ausgabedatei " << fileName << std::endl;
		return 1;
	}

	// Die VTK-Datei hat immer den gleichen Header. Dann schreiben wir die Koordinaten
	// der Punkte und anschließend die Definition der Line.
	output << "# vtk DataFile Version 3.0" << std::endl;
	output << "Linie aus PDEPack" << std::endl;
	output << "ASCII" << std::endl;
	output << "DATASET POLYDATA" << std::endl;
	output << "POINTS " << n << " float" << std::endl;
	// Die Koordinaten in einer Schleife
	// Wir tauschen die vorzeichen der y-Koordinaten, damit es zu den Bildern
	// im Text passt.
	for (i = 0; i < n; i++)
		output << std::fixed << x[i] << " " << -y[i] << " " << z << std::endl;
	output << "LINES 1 " << n + 1 << " " << std::endl;
	output << n << " ";
	for (i = 0; i < n; i++)
		output << i << " ";
	output << std::endl;
	output.close();
	return EXIT_SUCCESS;
}

short storeMatrix(string fileName,
	TNT::Matrix<double> values)
{
	short n(values.num_rows()), nSquare(n*n);
	double h(1.0 / static_cast<double>(n - 1));
	ofstream output;

	// open ostream
	output.open(fileName.c_str(), std::ios::out);
	if (output.fail()) {
		std::cerr << "Fehler beim Öffnen der Ausgabedatei " << fileName << std::endl;
		return 1;
	}

	// Die VTK-Datei hat immer den gleichen Header. Dann schreiben wir die Koordinaten
	// der Punkte und anschließend die Definition der Line.
	output << "# vtk DataFile Version 3.0" << std::endl;
	output << "Lösung einer Poisson-Gleichung" << std::endl;
	output << "ASCII" << std::endl;
	output << "DATASET STRUCTURED_POINTS" << std::endl;
	output << "DIMENSIONS " << n << " " << n << " 1" << std::endl;
	output << "ORIGIN 0 0 0" << std::endl;
	output << "SPACING " << h << " " << h << " 1" << std::endl << std::endl;
	output << "POINT_DATA " << nSquare << std::endl;
	output << "SCALARS Skalare double" << std::endl;
	output << "LOOKUP_TABLE default" << std::endl;
	// Die Koordinaten in einer Schleife
	// Die Daten werden zeilenweise ausgebeben
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			output << std::fixed << values[i][j] << " ";
	    output << std::endl;
    }
	output.close();
	return EXIT_SUCCESS;
}

