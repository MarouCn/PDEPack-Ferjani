#ifndef ITERATIVE_H
#define ITERATIVE_H

#include <tnt.h>

//! Rein virtuelle Basisklasse für iterative Algorithmen zur Lösung von linearen Gleichungssystemen
/*!
 * Die Klasse implementiert ein Default-Abbruchkriterium, das überschrieben werden kann. 
 * Die abgeleiteten Klassen implementieren die Funktion solve, die für eine übergebene
 * rechte Seite eine Lösung des linearen Gleichungssystems annähert.
 */
class Iterative
{
public:
    
        //! Defaultkonstruktor
	    /*!
	     *  n wird 0 gesetzt.
		 *  Abbruchgenauigkeit ist Eps aus utils.h
		 *
		 *  Maximale Anzahl der Iterationen ist 0
		 */
	    Iterative();

        //! Konstruktor mit Koeffizientenmatrix und einer rechten Seite
        /*!
		 *  Die Diagonallemente werden überprüft. Falls ein Diagonalelement zu klein ist gibt es eine
		 *  Ausgabe auf std::cerr. Dann sollte die Iteration nicht
		 *  durchgeführt werden.
		 *
		 * Abbruchgenauigkeit ist Eps aus utils.h
		 *
		 *  Maximale Anzahl der Iterationen ist 10
		 *
		 * \param rhs rechte Seite des Gleichungssystems
         */
		Iterative(TNT::Vector<double> rhs);

		//! Konstruktor mit Koeffizientenmatrix, einer rechten Seite, Abbruchgenauigkeit und maximaler Anzahl der Iterationen
		/*!
		 *  Die Diagonallemente werden überprüft. Falls ein Diagonalelement zu klein ist gibt es eine
		 *  Ausgabe auf std::cerr. Dann sollte die Iteration nicht
		 *  durchgeführt werden.
		 *
		 * \param rhs rechte Seite des Gleichungssystems
		 * \param epsilon Abbruchgenauigkeit
		 * \param maxI Maximale Anzahl der Iterationen
		 */
		Iterative(TNT::Vector<double> rhs, double epsilon, unsigned int maxI);

        //! Lösung des linearen Gleichungssystems
	    /*!
		 * Generische Durchführung des iterativen Verfahrens. Verwendet
		 * die abstrakte Funktion Iterative::iteration().
		 *
		 * Überprüft, ob die maximale Anzahl der Iterationen bereits erreicht
		 * wurde. Das Abbruchkriterium muss in der Funktion implementiert werden,
		 * die den Iterationsschritt implementiert.
		 *
		 * Diese Funktion ist generisch und kann nicht überschrieben werden!
		 *
		 * \return Aktuelle Näherung wird zurückgegeben.
	     */
		TNT::Vector<double> solve(TNT::Vector<double> start);

		//! Durchführen eines Iterationsschritts
		/*!
		* Funktion, in der ein Schritt des Algorithmus durchgeführt wird.
		* In dieser Funktion wird das Abbruchkriterium durchgeführt
		* und, falls es erfüllt ist, muss die Funktion die Variable
		* Iterative::shouldStop auf true gesetzt werden. 
		*/
		virtual void iterate() = 0;

        //! Ausgabe der Inhalt der Instanz auf der Konsole
        virtual void print() const;

		//! Abfragen einer Lösung
		TNT::Vector<double> getSolution();

		//! Setzen der rechten Seite
		void setRhs(TNT::Vector<double> rhs);

		//! Setzen der Abbruchgenauigkeit
		void setEpsilon(double e);

		//! Setzen der maximalen Anzahl der Iterationen
		void maxNumberOfIterations(unsigned m);

protected:
		//! Dimension des Gleichungssystems
		unsigned int n;
		//! Vektor mit der rechten Seite des linearen Gleichungssystems
		TNT::Vector<double> b;
		//! Feld für die Durchführung der Iteration
		/*!
		 * Das Feld enthält nach der Iteration die berechnete Lösung und
		 * kann mit Iterative::getSolution() abgefragt werden.
		 */
		TNT::Vector<double> x;
		//! Abbruchgenauigkeit
		double localEpsilon;
		//! Maximale Anzahl der Iterationen als "Notbremse"
		unsigned int maxIter;
		//! Anzahl der Iterationen
		/*!
		 * Nach dem Aufruf von Iterative::solve kann mit dieser Variable die benötigte
		 * Anzahl der Iterationen, die durchgeführt wurden, ausgegeben werden.
		 * Vorher ist dieser Wert 0.
		 */
		unsigned int noI;
		//! Falls diese Variable true ist ist das Abbruchkriterium erfüllt
		bool shouldStop;
		//! Zeigt an, dass das Verfahren beendet ist
		/*!
		 * Das Verfahren wurde beendet, entweder auf Grund des Abbruchkriteriums
		 * oder auf Grund des Erreichens der maximalen Anzahl der Iterationen.
		 */
		bool finished;
		//! Funktion,die die Iteration nach Erfüllung des Abbruchkriteriums anhält
		void stop();
		//! Zurücksetzen der Parameter, um eine erneute Berechnung zu ermöglichen
		void reset();
};

// - Inline section ---------------------------------------

// Setzen der Abbruchgenauigkeit, wenn nicht der Default verwendet werden soll
inline void Iterative::setEpsilon(double e) { localEpsilon = e; }
// Setzen der maximalen Anzahl der Iterationen
inline void Iterative::maxNumberOfIterations(unsigned int m) { maxIter = m;  }
// Wechseln der rechten Seite
inline void Iterative::setRhs(TNT::Vector<double> rhs) { b = rhs; }
// Anhalten der Iteration
inline void Iterative::stop() { shouldStop = true; }
// Rücksetzebn der Iteration
inline void Iterative::reset() { shouldStop = false; finished = false; noI = 0; }
// ALösung abfragen
inline TNT::Vector<double> Iterative::getSolution() { return x; }

#endif
