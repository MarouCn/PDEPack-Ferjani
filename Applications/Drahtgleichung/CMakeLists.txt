# --------------------------------------------------------
#  Projekt für LU-Zerlegung von tridiagonalen Matrizen
# --------------------------------------------------------
# LU-Zerlegung
add_executable (Drahtgleichung EXCLUDE_FROM_ALL drahtgleichung.cpp)
TARGET_LINK_LIBRARIES(Drahtgleichung pdepack)
