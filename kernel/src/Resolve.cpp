#include "Resolve.h"
/*!
 * \brief Resolve::Resolve
 * In dieser Klasse sind die 3 Algorithmus implementiert(Jacobi , Gauss Seidel , Successive Over-Relaxation )
 * Die VTK Datein fuer Paraview werden auch mit der Hilfe dieser Klasse erzeugt
 * Benchmark Datein werden dazu erstellt
 */
Resolve::Resolve()
{
	this->matrix = TNT::Matrix<double>();
	this->dimension = 0;
	this->temperature = 0;
	this->f = -2.0;
	this->h = 0;
	this->useOpenMP = false;
}


Resolve::Resolve(int temperature, int dimension, int timeLimit, bool useOpenMP)
{
        try
        {
                this->matrix = TNT::Matrix<double>(dimension, dimension);
        }
        //Speicher kann nicht allokiert werden, wahrscheinlich ist die Matrix zu groﬂ
        catch (std::bad_alloc)
        {
                std::cout << "Nicht genug Arbeitsspeicher zur Verfuegung um eine Matrix dieser Dimension zu erstellen. Es wird stattdessen die Dimension 250 verwendet." << std::endl;
                this->matrix = TNT::Matrix<double>(250, 250);
                dimension = 250;
        }
        this->temperature = temperature;
        this->dimension = dimension;
        this->f = -2.0;
        this->h = (double)1 / dimension;
        this->timeLimit = timeLimit;
        this->useOpenMP = useOpenMP;
}




void Resolve::calculateSOR()
{
	algo = "SOR";
	//optimalen Relaxationsparameter berechnen
	double woptimal = 2 / (1 + sqrt(1 - pow(cos(M_PI / (dimension + 1)), 2)));
	bool done = false;

	//Temperatur in die Matrix einf¸gen
	for (int i = 0; i < dimension; i++)
	{
		matrix[i][0] = temperature;
	}

	TNT::Matrix<double> matrixOld;
	try
	{
		matrixOld = *new TNT::Matrix<double>(matrix);
	}
	//Speicher kann nicht allokiert werden, wahrscheinlich ist die Matrix zu groﬂ
	catch (std::bad_alloc)
	{
		std::cout << "Nicht genug Arbeitsspeicher zur Verfuegung. Programm wird beendet." << std::endl;
		system("pause");
	}

	std::chrono::time_point<std::chrono::high_resolution_clock> startTime = getWalltime();
	double consumedTime = 0.0;

	//Berechnung
	for (count = 0; !done && (timeLimit == 0 || consumedTime < timeLimit); count++)
	{
		#pragma omp parallel for if (useOpenMP)
		for (int i = 1; i < dimension - 1; i++)
		{
			#pragma omp parallel for if (useOpenMP)
			for (int j = 1; j < dimension - 1; j++)
			{
				matrix[i][j] = (1 - woptimal) * matrix[i][j] + (woptimal / 4) *(matrix[i - 1][j] + matrix[i + 1][j] + matrix[i][j - 1] + matrix[i][j + 1] + (-(h*h)*f));
			}
		}
		done = compareMatrixNorm(matrixOld, matrix, Eps);
		matrixOld = matrix;
		consumedTime = elapsedSeconds(startTime, getWalltime()).count();
	}
	//Berechnungsdauer
	duration = elapsedSeconds(startTime, getWalltime()).count();
}

void Resolve::calculateGS()
{
	algo = "Gauss-Seidel";
	bool done = false;

	//Temperatur einf¸gen
	for (int i = 0; i < dimension; i++)
	{
		matrix[i][0] = temperature;
	}

	TNT::Matrix<double> matrixOld;
	try
	{
		matrixOld = *new TNT::Matrix<double>(matrix);
	}
	//Speicher kann nicht allokiert werden, wahrscheinlich ist die Matrix zu groﬂ
	catch (std::bad_alloc)
	{
		std::cout << "Nicht genug Arbeitsspeicher zur Verfuegung. Programm wird beendet." << std::endl;
		system("pause");
	}

	std::chrono::time_point<std::chrono::high_resolution_clock> startTime = getWalltime();
	double consumedTime = 0.0;

	//Berechnung
	for (count = 0; !done && (timeLimit == 0 || consumedTime < timeLimit); count++)
	{
			#pragma omp parallel for if (useOpenMP)
			for (int i = 1; i < dimension - 1; i++)
			{
				#pragma omp parallel for if (useOpenMP)
				for (int j = 1; j < dimension - 1; j++)
				{
					matrix[i][j] = (1.0 / 4.0)*(matrix[i - 1][j] + matrix[i + 1][j] + matrix[i][j - 1] + matrix[i][j + 1] + (-(h*h)*f));
				}
			}
			done = compareMatrixNorm(matrixOld, matrix, Eps);
			matrixOld = matrix;
		    consumedTime = elapsedSeconds(startTime, getWalltime()).count();
	}
	//Berechnungsdauer
	duration = elapsedSeconds(startTime, getWalltime()).count();
}

void Resolve::calculateJ()
{
	algo = "Jacobi";
	bool done = false;

	//Temperatur einf¸gen
	for (int i = 0; i < dimension; i++)
	{
		matrix[i][0] = temperature;
	}

	TNT::Matrix<double> matrixOld;
	try
	{
		matrixOld = *new TNT::Matrix<double>(matrix);
	}
	//Speicher kann nicht allokiert werden, wahrscheinlich ist die Matrix zu groﬂ
	catch (std::bad_alloc)
	{
		std::cout << "Nicht genug Arbeitsspeicher zur Verfuegung. Programm wird beendet." << std::endl;
		system("pause");
	}

	std::chrono::time_point<std::chrono::high_resolution_clock> startTime = getWalltime();
	double consumedTime = 0.0;

	//Berechnung
	for (count = 0; !done && (timeLimit == 0 || consumedTime < timeLimit); count++)
	{
		#pragma omp parallel for if(useOpenMP)
		for (int i = 1; i < dimension - 1; i++)
		{
			#pragma omp parallel for if(useOpenMP)
			for (int j = 1; j < dimension - 1; j++)
			{
				matrix[i][j] = (1.0 / 4.0) * (matrixOld[i - 1][j] + matrixOld[i + 1][j] + matrixOld[i][j - 1] + matrixOld[i][j + 1] + (-(h*h)*f));
			}
		}
		done = compareMatrixNorm(matrixOld, matrix, Eps);
		matrixOld = matrix;
		consumedTime = elapsedSeconds(startTime, getWalltime()).count();
	}
	std::cout << useOpenMP;
	//Berechnungsdauer
	duration = elapsedSeconds(startTime, getWalltime()).count();
}


bool Resolve::exportToVTK(std::string fileName) const
{

	std::ofstream outStream;
	outStream.open(fileName);
	if (outStream.fail()) return false;

	//Header schreiben
	outStream << "# vtk DataFile Version 2.0" << std::endl;
	outStream << "Temperaturverlauf in einer K¸hlrippe mit konstanter W‰rmequelle von " << temperature << " Grad Celsius und Matrix-Dimension " << dimension << "x" << dimension << std::endl;
	outStream << "ASCII" << std::endl;
	outStream << "DATASET STRUCTURED_POINTS" << std::endl;
	outStream << "DIMENSIONS " << dimension << " " << dimension << " 1" << std::endl;
	outStream << "ORIGIN 0 0 0" << std::endl;
	outStream << "SPACING " << h << " " << h << " 1" << std::endl;
	outStream << "" << std::endl;
	outStream << "POINT_DATA " << dimension * dimension << std::endl;
	outStream << "SCALARS Skalare double" << std::endl;
	outStream << "LOOKUP_TABLE default" << std::endl;

	//Matrixwerte schreiben
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			outStream << matrix[i][j] << " ";
		}
		
		outStream << std::endl;
	}
	

	outStream.close();
	return true;
}

bool Resolve::appendToLog(std::string fileName) const
{

	std::ofstream logFile;
	logFile.open(fileName, std::ios_base::app);
	if (logFile.fail()) return false;
	//Daten schreiben
        logFile << "Temperatur: " << temperature << ", Dimension: " << dimension << ", f(xi,xj): " << f << ", Algorithmus: " << algo  << " useOpenMP: " << useOpenMP << ":" << std::endl << "Dauer der Berechnung: " << duration << "ms, Iterationen: " << count << std::endl << std::endl;

	logFile.close();
	return true;
}

int Resolve::getIterations()
{
	return count;
}

double Resolve::getDuration()
{
	return duration;
}

TNT::Matrix<double> Resolve::getMatrix()
{
	return matrix;
}
