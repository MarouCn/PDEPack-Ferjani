#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "utils.h"
#include "Wire1D.h"  

Wire1D::Wire1D() : l(0.0), n(0), delta(0.0), a(0.0), b(0.0) {}

Wire1D::Wire1D(int np) : l(1.0), n(np), delta(l/static_cast<double>(n+1)),
        a(0.0), b(0.0), 
		tri(n, 2.0, -1.0, -1.0), 
		lu(tri) {}


Wire1D::Wire1D(double len, unsigned int np, double left, double right) 
	 : l(len), n(np), delta(l / static_cast<double>(n + 1)),
               a(left), b(right),
               tri(n, 2.0, -1.0, -1.0),
               lu(tri) 
{}

TNT::Vector<double> Wire1D::solve(TNT::Vector<double> forces) const
{
	TNT::Vector<double> rhs(forces);

	for (unsigned int i = 0; i < n; i++)
		rhs[i] *= delta*delta;
	rhs[0] += a;
	rhs[n-1] += b;

	TNT::Vector<double> sol(lu.solve(rhs));

	// Ausgabe des Kabels, inkl. der Randpunkte
	TNT::Vector<double> y(n + 2);
	y[0] = a;
	for (unsigned int i = 1; i < n + 1; i++)
		y[i] = sol[i-1];
	y[n + 1] = b;
	return y;
}

TNT::Vector<double> Wire1D::getXCoord() const
{
	TNT::Vector<double> x(n + 2);
	x[0] = 0.0;
	for (unsigned int i = 1; i < n + 1; i++)
		x[i] = static_cast<double>(i)* delta;
	x[n + 1] = l;
	return x;
}