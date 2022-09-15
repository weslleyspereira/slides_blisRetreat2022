cmake_minimum_required(VERSION 3.11)
project( testsBLISRetreat C )

find_package( BLAS REQUIRED )

# add the example eigen
add_executable( testBLIS testBLIS.c )
target_link_libraries( testBLIS PRIVATE ${BLAS_LIBRARIES} -lm )