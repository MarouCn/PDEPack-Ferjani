# --------------------------------------------------------
#  Projekt für LU-Zerlegung von tridiagonalen Matrizen
# --------------------------------------------------------
# LU-Zerlegung
add_executable (luTest EXCLUDE_FROM_ALL luTest.cpp)
TARGET_LINK_LIBRARIES(luTest pdepack)

# LU-Zerlegung
add_executable (triTest EXCLUDE_FROM_ALL triTest.cpp)
TARGET_LINK_LIBRARIES(triTest pdepack)
