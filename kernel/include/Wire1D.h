#ifndef WIRE1D_H
#define WIRE1D_H

#include <tnt.h>

#include "TriDiagScalar.h"
#include "TriLUScalar.h"

//! Klasse für die numerische Lösung der eindimensionalen Kabelgleichung
class Wire1D 
{
public:
        //! Defaultkonstruktor
	    /*!
	     *  n w=0
		 *  delta=0.0
		 *  a = b = 0
		 */
	    Wire1D();

		//! Konstruktor mit homogenen Randbedingungen und Länge 1
		/*!
		 *
		 *  Die Randwerte sind 0.0, die Länge ist 1.
		 *
		 * \param np Anzahl der Punkte in der Diskretisierung
		 */
		Wire1D(int np);
        //! Konstruktor mit allen Angaben zur Drahtgleichung
        /*!
         * \param len Länge
		 * \param np Anzahl der Punkte in der Diskretisierung
		 * \param left Randwert bei x=0
		 * \param right Randwert bei x=l
         */
		Wire1D(double len, unsigned int np, double left=0.0, double right=0.0);
	
		//! Berechnung einer Näherung für den aktuellen Kraftvektor
		/*!
		 * \return Aktuelle Näherung wird zurückgegeben, incl. Randwerte
		 */
		TNT::Vector<double> solve(TNT::Vector<double> forces) const;

		//! Ausgabe der x-Koordinaten er Diskretisierung des Kabels
		/*!
		 * Ausgegeben werden n+2 Werte, die innere Punkte, an
		 * denen die Lösung der partiellen Differenzialgleichung
		 * berechnet wird und die zwei Randpunkte x=0.0
		 * und x=l.
		 *
		 * Die dazu gehörigen y-Koordinaten erhalten wir als Ausgabe
		 * der Funktion Wire1D::solve;
		 */
		TNT::Vector<double> getXCoord() const;

private:

	//! Länge des Kabels
	/*!
	 * Wir gehen davon aus, dass das "linke" Ende des Kabels im Ursprung liegt,
	 * und das Kabel in positiver x-Richtung verläuft.
	 */
	double l;
	//! Anzahl der Punkte in der Diskretisierung
	/*!
	 * Zu den n Werten kommen noch die beiden Randwerte dazu.
	 */
	unsigned int n;
	//! Schrweitweite in der Diskretisierung
	double delta;
	//! Randwert für x=0 
	double a;
	//! Randwert für x=l
	double b;
	//! Tridiagonalmatrix für die numerische Lösung
	TriDiagScalar tri;
	//! Instanz einer LU-Zerlegung für tridiagonale Matrizen
	TriLUScalar lu;
};

#endif
