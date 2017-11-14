#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include <tnt.h>
#include "Iterative.h"

//! Einzelschritt-Verfahren zur Lösung eines linearen Gleichungssystems
/*!
 * Diese Klasse realisiert die Lösung eines linearen Gleichungssytems
 * mit Hilfe des Einzelschritt- oder Gauss-Seidel-Verfahrens.
 *
 * Es wird davon ausgegangen, dass die Koeffizientenmatrix eine voll
 * besetzte Matrix ist.
 */
class GaussSeidel : public Iterative
{
public:
        //! Defaultkonstruktor
	    /*!
	     *  n wird 0 gesetzt.
		 *  Abbruchgenauigkeit ist Eps aus utils.h
		 *
		 *  Maximale Anzahl der Iterationen ist 0
		 */
	    GaussSeidel();

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
		GaussSeidel(TNT::Matrix<double> matrix, TNT::Vector<double> rhs);

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
		GaussSeidel(TNT::Matrix<double> matrix, TNT::Vector<double> rhs, double epsilon, unsigned int maxI);

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
	//! Zeiger auf eine Instanz einer  Matrix
	TNT::Matrix<double> A;
	//! Vektor für die vorhergehende Lösung in der Iteration
	/*!
	 * Dieser Vektor wird ausschließlich für das Abbruchkriterium benötigt!
	 */
	TNT::Vector<double> xalt;
};

#endif
