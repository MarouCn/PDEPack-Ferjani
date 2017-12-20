# --------------------------------------------------------
#  Loesung Programm fuer die Temperature Distribution in einer Kühlrippe
# --------------------------------------------------------
add_executable (kuehlrippeTest EXCLUDE_FROM_ALL kuehlrippeTest.cpp)
TARGET_LINK_LIBRARIES(kuehlrippeTest pdepack)