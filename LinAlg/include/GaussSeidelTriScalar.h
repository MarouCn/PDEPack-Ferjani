#ifndef GAUSSSEIDELTRISCALAR_H
#define GAUSSSEIDELTRISCALAR_H

#include <tnt.h>
#include "Iterative.h"
#include "TriDiagScalar.h"

//! Einzelschritt-Verfahren zur Lösung eines linearen Gleichungssystems
/*!
 * Diese Klasse realisiert die Lösung eines linearen Gleichungssytems
 * mit Hilfe des Einzelschritt- oder GaussSeidel-Verfahrens.
 *
 * Die Koeffizientenmatrix ist eine tridiagonale Matrix mit
 * konstanten Werten auf den Diagonalen, realisiert durch die Klasse
 * TriDiagScalar.
 */
class GaussSeidelTriScalar : public Iterative
{
public:
    
        //! Defaultkonstruktor
	    /*!
	     *  n wird 0 gesetzt.
		 *  Abbruchgenauigkeit ist Eps aus utils.h
		 *
		 *  Maximale Anzahl der Iterationen ist 0
		 */
	    GaussSeidelTriScalar();

		//! Konstruktor mit den Diagonalelementen
		/*!
		*  Falls ein Diagonalelement zu klein ist gibt es eine
		*  Ausgabe auf std::cerr. Dann sollte die Iteration nicht
		*  durchgeführt werden.
		*
		* Abbruchgenauigkeit ist Eps aus utils.h
		*
		*  Maximale Anzahl der Iterationen ist 10
		*
		* \param d Diagonale
		* \param up Obere Nebendiagonale
		* \param low Untere Nebendiagonale
		* \param rhs rechte Seite des Gleichungssystems
		*/
		GaussSeidelTriScalar(double d, double up, double low, TNT::Vector<double> rhs);

        //! Konstruktor mit einer Koeffizientenmatrix
        /*!
		 *  Falls ein Diagonalelement zu klein ist gibt es eine
		 *  Ausgabe auf std::cerr. Dann sollte die Iteration nicht
		 *  durchgeführt werden.
		 *
		 * Abbruchgenauigkeit ist Eps aus utils.h
		 *
		 *  Maximale Anzahl der Iterationen ist 10
		 *
         * \param matrix Koeffizientenmatrix
		 * \param rhs rechte Seite des Gleichungssystems
         */
		GaussSeidelTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs);

		//! Konstruktor mit Koeffizientenmatrix, einer rechten Seite, Abbruchgenauigkeit und maximaler Anzahl der Iterationen
		/*!
		*  Die Diagonallemente werden überprüft. Falls ein Diagonalelement zu klein ist gibt es eine
		*  Ausgabe auf std::cerr. Dann sollte die Iteration nicht
		*  durchgeführt werden.
		*
		* \param matrix Koeffizientenmatrix
		* \param rhs rechte Seite des Gleichungssystems
		* \param epsilon Abbruchgenauigkeit
		* \param maxI Maximale Anzahl der Iterationen
		*/
		GaussSeidelTriScalar(TriDiagScalar matrix, TNT::Vector<double> rhs, double epsilon, unsigned int maxI);

		//! Durchführen eines Iterationsschritts
		/*!
		* Funktion, in der ein Schritt des Algorithmus durchgeführt wird.
		* In dieser Funktion wird das Abbruchkriterium durchgeführt
		* und, falls es erfüllt ist, muss die Funktion die Variable
		* Iterative::shouldStop auf true gesetzt werden.
		*
		* \return Aktuelle Näherung wird zurückgegeben.
		*/
		virtual void iterate();

        //! Ausgabe der Inhalt der Instanz auf der Konsole
        void print() const;
private:
	    //! Diagonalelement
		double diag;
		//! Oberes Diagonalelement
		double upper;
		//! Unteres Diagonalelement
		double lower;
	    //! Vektor für die vorhergehende Lösung in der Iteration
		/*!
		 * Wir für das Einzelschrittverfahren nur für das Abbruchkriterium benötigt!
		 */
	    TNT::Vector<double> xalt;

};

// - Inline section ---------------------------------------


#endif
