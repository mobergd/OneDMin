project(OneDMin Fortran)
cmake_minimum_required(VERSION 2.8)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(ljlibs
    ${PROJECT_SOURCE_DIR}/src/mass.f
    ${PROJECT_SOURCE_DIR}/src/spin.f)

add_library(exp6libs
    ${PROJECT_SOURCE_DIR}/src/exp6.f)

add_library(molprolibs
    ${PROJECT_SOURCE_DIR}/src/molpro.f)

message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_executable(onedmin-exp6.x
    ${PROJECT_SOURCE_DIR}/src/onedmin.f)

target_link_libraries(onedmin-exp6.x 
    ljlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} exp6libs)

add_executable(onedmin-molpro.x
    ${PROJECT_SOURCE_DIR}/src/onedmin.f)

target_link_libraries(onedmin-molpro.x 
    ljlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} molprolibs)

install(TARGETS onedmin-exp6.x DESTINATION bin)
install(TARGETS onedmin-molpro.x DESTINATION bin)
