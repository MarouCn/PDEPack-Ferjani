# --------------------------------------------------------
#  Loesung Programm fuer die Temperature Distribution in einer K�hlrippe
# --------------------------------------------------------
add_executable (kuehlrippeTest EXCLUDE_FROM_ALL kuehlrippeTest.cpp)
TARGET_LINK_LIBRARIES(kuehlrippeTest pdepack)