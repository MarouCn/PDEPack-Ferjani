#include <cmath>
#include "utils.h"


// Funktion für die Besetzung des i-ten kanonischen Einheitsvektors im R^n
TNT::Vector<double> unitVector(int i, int n)
{
	TNT::Vector<double> e(n);
	e[i] = 1.0;
	return e;
}

// Vergleich eines TNT Vektors mit einem Zielvektor mit Hilfe der euklidischen Norm
// Ist die relative Abweichung kleiner als die Konstante Eps, dann wird true zurückgegeben
bool compareVector(TNT::Vector<double> value, TNT::Vector <double> target, double epsilon)
{
	TNT::Vector<double> diff(value - target);
	return (sqrt(TNT::dot_prod(diff, diff) / TNT::dot_prod(target, target)) < epsilon);
}

// Vergleich eines TNT Vektors mit einem Zielvektor mit Hilfe der Maximums-Norm
// Ist die relative Abweichung kleiner als die Konstante Eps, dann wird true zurückgegeben
bool compareVectorMaximumNorm(TNT::Vector<double> value, TNT::Vector<double> target, double epsilon)
{
	unsigned int n(value.size());
	for (unsigned int i = 0; i<n; i++)
	{
		if (fabs((value[i] - target[i]) / target[i])>epsilon)
		{
			return false;
		}
	}
	return true;
}

// Berechnung der Maximumsnorm eines Vektors
double normMaximum(TNT::Vector<double> v)
{
	unsigned int n(v.size());
	double max(0.0);
	for (unsigned int i = 0; i < n; i++)
	{
		if (fabs(v[i]) > max) {
			max = fabs(v[i]);
		}
	}
	return max;
}