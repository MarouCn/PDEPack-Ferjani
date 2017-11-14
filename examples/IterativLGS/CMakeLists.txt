# --------------------------------------------------------
#  Iterative Verfahren für LGS
# --------------------------------------------------------
# Jacobi - Gesamtschrittverfahren
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/matrix1.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/matrix2.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/matrix3.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_executable (jacobiTest EXCLUDE_FROM_ALL jacobiTest.cpp matrix1.txt)
add_executable (gaussSeidelTest EXCLUDE_FROM_ALL gaussSeidelTest.cpp matrix1.txt)
add_executable (triJacobiTest EXCLUDE_FROM_ALL triJacobiTest.cpp matrix3.txt)
add_executable (triGaussSeidelTest EXCLUDE_FROM_ALL triGaussSeidelTest.cpp matrix3.txt)

TARGET_LINK_LIBRARIES(jacobiTest pdepack)
TARGET_LINK_LIBRARIES(gaussSeidelTest pdepack)
TARGET_LINK_LIBRARIES(triJacobiTest pdepack)
TARGET_LINK_LIBRARIES(triGaussSeidelTest pdepack)

