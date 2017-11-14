# ------------------------------------------------------------
#  CMakeLists.txt für die Fallstudie PDEPack
# ------------------------------------------------------------
cmake_minimum_required (VERSION 3.6)
project (PDEPack)

# Generelle Optionen, File wird auch in Examples eingelesen!
include (CMakeOptions.txt)
# Einstellungen für Doxygen
include (CMakeDoxygen.txt) 

# Alle Quellen und Header zur Bibliothek hinzufügen
# Die Header müssen explizit genannt werden damit sie zum Microsoft Visual Studio
# Projekt hinzugefügt werden!
FILE(GLOB SRCS_CPP "kernel/src/*.cpp" "LinAlg/src/*.cpp" "Vis/src/*.cpp")
# Die folgende Zeile erzeugt eine shared library, wenn das Kommentarzeichen entfernt wird!
FILE(GLOB SRCS_H "kernel/include/*.h" "LinAlg/include/*.h" "Vis/include/*.h")

# Sicher stellen dass die Header gefunden werden
include_directories (${PDEPack_SOURCE_DIR}/kernel/include ${PDEPack_SOURCE_DIR}/LinAlg/include ${PDEPack_SOURCE_DIR}/Vis/include)

# Die Bibliotheken erzeugen
# Wir erzeugen sowohl eine statische als auch eine dynamische Bibliothek
add_library (pdepack STATIC ${SRCS_CPP} ${SRCS_H})
# Wenn man die folgende Zeile dekommentiert erhält man eine shared library
#add_library (pdepack_shared SHARED ${SRCS_CPP} ${SRCS_H})

# Beispiele in das Projekt einbeziehen oder nicht
option(EXAMPLES "Beispielprojekte erstellen" OFF)
if(EXAMPLES)
     add_subdirectory(examples)
endif(EXAMPLES)

# Anwendungen im Verzeichnis Applications in das Projekt einbeziehen oder nicht
option(APPS "Anwendungen erstellen" OFF)
if(APPS)
	add_subdirectory(Applications)
endif(APPS)
