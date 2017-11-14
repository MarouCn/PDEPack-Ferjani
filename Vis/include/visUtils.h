#ifndef VISUTIL
#define VISUTIL

#include <tnt.h>

//! \addtogroup visutil visutil - Nützliches, Konstanten, ...
/*! \addtogroup visutil
  *  @{
  */
  
//! Ausgabe einer Linie als VTK-Datei
/*!
 * Die VTK-Datei kann mit ParaView oder dem VTK-Framework
 * geöffnet und damit eine Visualisierung erstellt werden.
 *
 * Verwenden Sie die ParaView-Pipeline-Datei, die online
 * zur Verfügung steht. Oder öffnen Sie im Datei-Dialog
 * in Paraview die mit dieser Funktion erzeugte
 * vtk-Datei und fügen Sie eine Farbtabelle hinzu.
 *
 * \param fileName Dateiname
 * \param x Vektor mit den x-Koordinaten der Diskretisierung
 * \param y Vektor mit den y-Koordinaten der Diskretisierung
 * \return Fehlercodes
 *
 * 0 = success
 * 1 = Fehler beim Öffnen der Ausgabedatei
 */
short storeLine(string fileName, 
	            TNT::Vector<double> x, 
				TNT::Vector<double> y);

//! Ausgabe einer zweidimensionalen Datenstruktur als VTK-Datei
/*!
* Die VTK-Datei kann mit ParaView oder dem VTK-Framework
* geöffnet und damit eine Visualisierung erstellt werden.
*
* Verwenden Sie die ParaView-Pipeline-Datei, die online
* zur Verfügung steht. Oder öffnen Sie im Datei-Dialog
* in Paraview die mit dieser Funktion erzeugte
* vtk-Datei und fügen Sie eine Farbtabelle hinzu.
*
* \param fileName Dateiname
* \param values Matrix mit den Werten an den n*n Punkten, also einschließlich
*               der Randwerte! Die Werte werden zeilenweise erwartet. Bei 
*               n=8 werden 8 Zeilen mit jeweils 8 Einträgen in die Datei
*               beschrieben.
*               Es müssen alle Werte, auch die Randwerte, übergeben werden.
*               Im stationären Fall mit 6 Punkten pro Zeile, für die die Näherung berechnet wird,
*               wird eine 8x8-Matrix erwartet.
* \return Fehlercodes
*
* 0 = success
* 1 = Fehler beim Öffnen der Ausgabedatei
*/
short storeMatrix(string fileName,
	              TNT::Matrix<double> values);

 /*! @}  end doxygen addtogroup visutil */	
#endif // VISUTIL
