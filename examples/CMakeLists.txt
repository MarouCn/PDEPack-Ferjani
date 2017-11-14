# CMake-File für die Hauptprogramme

# Name des Projekts: "PDEPack". 
#cmake_minimum_required (VERSION 3.4)
#project (PDEPack_Mains)

include (../CMakeOptions.txt)

# Jedes Beispiel liegt in einem Unterverzeichnis!
# Ein neues Beispiel muss mit add_subdirectory hinzugefügt werden!
add_subdirectory(TriLU)

add_subdirectory(IterativLGS)
