#ifndef TRILUSCALAR_H
#define TRILUSCALAR_H


#include "TriLU.h"
#include "TriDiagScalar.h"

//! LU-Zerlegung und Auflösung von linearen Gleichungssystemen mit einer tridiagonalen Matrix, die nur konstante Werte in den Diagonalen besitzt
/*!
 * Wir implementieren ausschließlich den Fall, dass die Berechnung der LU-Zerlegung
 * ohne Zeilenvertauschung durchführbar ist.
 * Dies ist in sehr vielen Anwendungen garantiert.
 *
 * Wir verwenden die gleiche Notation wie in der Musterlösung der Übungsaufgabe. 
 *
 * Die Diagonale des Faktors L ist durch 1 gegeben.
 */
class TriLUScalar : public TriLU
{
public:
    
        //! Defaultkonstruktor
	    /*!
	     *  n wird 0 gesetzt.
		 */
	    TriLUScalar();

		//! Konstruktor mit Dimension der linearen Gleichungssystems und konstanten Werten für die Diagonalen
		/*!
		 * \param lda Dimensiona des Gleichungssystems
		 * \param diag Konstantes Element auf der Diagonalen der Koffizientenmatrix
		 * \param upper Konstantes Element auf der oberen Nebendiagonale der Koeffizientenmatrix
		 * \param lower Konstantes Element auf der unteren Nebendiagonale der Koeffizientenmatrix
		 */
		TriLUScalar(unsigned int lda, double diag, double upper, double lower);

        //! Konstruktaor mit einer Koeffizientenmatrix
        /*!
		 * Analog zur TNT::LU Zerlegung wird die LU-Zerlegung direkt 
		 * im Konstruktor berechnet.
		 *
         * \param matrix Instanz einer tridiagonalen Matrix mit konstanten Werten auf der Diagonale
         */
		TriLUScalar(TriDiagScalar matrix);

        //! Ausgabe der Inhalt der Instanz auf der Konsole
        void print() const;


protected:

		//! Berechnung der LU-Zerlegung ohne Zeilenvertauschung
		/*!
		* @return 0 Elimination erfolgreich
		*         1 Elimination abgebrochen, Diagonalelement zu klein
		*/
		virtual short decompose();

	   //! Eine Instanz einer tridiagonalen Matrix
	   TriDiagScalar M;
};
#endif
