#ifndef TRIDIAGARRAY_H
#define TRIDIAGARRAY_H

#include <tnt.h>

//! Repräsentation einer tridiagonalen Matrix mit Hilfe von drei TNT-Vektoren
class TriDiagArray 
{

public:
    
        //! Defaultkonstruktor
	    /*!
	     *  n wird 0 gesetzt, die Felder für die Koeffizientenmatrix
	     *  ebenfalls!
	     */
        TriDiagArray();

		//! Konstruktor mit der Dimension der Matrix
		/*!
		 *  Nur die Dimension wird besetzt, die Felder werden nicht angelegt!
		 */
		TriDiagArray(unsigned int lda);
        //! Konstruktor mit Arrays für die Diagonalelemente
        /*!
		 * Für eine nxn tridiagonale Matrix hat die Diagonale
		 * lda Elemente, die beiden Nebendiagonalen lda-1 Elemente!
		 *
         * \param diag Die lda Elemente auf der Diagonale
	     * \param upperDiag Die lda-1 Elemente auf der oberen Nebendiagonale
	     * \param lowerDiag Die lda-1 Elemente auf der unteren Nebendiagonale
	     *
	     * Wichtig: es wird nicht überprüft, ob die übergebenen
	     * Felder die korrekte Dimensionen haben!
         */
		TriDiagArray(TNT::Vector<double> diag, TNT::Vector<double> upperDiag, TNT::Vector<double> lowerDiag);

        //! Konstruktor für eine Matrix mit konstanten Diagonal-Elementen
        /*!
         * \param lda Anzahl Zeilen und Spalten
         * \param a Der konstante Wert für die Diagonale
	     * \param b Der konstante Wert für die obere Nebendiagonale
	     * \param c Der konstante Wert für die untere Nebendiagonale
         */
        TriDiagArray(unsigned int lda, double a, double b, double c);

        //! Multiplizieren der Matrix mit einem übergebenen Vektor
        /*!
         * \param x Der Vektor für den A*x berechnet wird.
	     * \return Zeiger auf ein Feld der Länge lda mit dem Ergebnis der Multiplikation A*x
         */
		virtual TNT::Vector<double> ATimesX(TNT::Vector<double> x) const;

        //! Ausgabe der Inhalt der Instanz auf der Konsole
        virtual void print() const;

		//! Die Dimension des linearen Gleichungssystems abfragen
		/*!
		* \return Anzahl der Zeilen und Spalten der quadratischen Koeffizientenmatrix
		*/
		int getDimension() const;

		//! Abfragen der Diagonale
		virtual TNT::Vector<double> getDiagonal() const;

		//! Abfragen der Elemente der unteren Nebendiagonale der Matrix
		virtual TNT::Vector<double> getUpperDiagonal() const;

		//! Abfragen des i-ten Elements der unteren Nebendiagonale
		virtual TNT::Vector<double> getLowerDiagonal() const;

protected:
	    //! Anzahl der Gleichungen und der Unbekannten
	    unsigned int n;

private:
        //! Vektor aus dem R^n mit den Diagonalelementen
	    TNT::Vector<double> a;
        //! Vektor aus dem R^(n-1) mit der unteren Nebendiagonale
	    TNT::Vector<double> c;
        //! Vektor aus dem R^(n-1) mit der oberen Nebendiagonale
	    TNT::Vector<double> b;
};

// - Inline section ---------------------------------------

// Abfragen der Dimension als inline-Funktion
inline int TriDiagArray::getDimension() const { return n; }
#endif
