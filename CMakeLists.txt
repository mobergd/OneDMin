project(Auto1DMin Fortran)
cmake_minimum_required(VERSION 2.8)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(ljlibs
    ${PROJECT_SOURCE_DIR}/src/molpro.f
    ${PROJECT_SOURCE_DIR}/src/mass.f
#    ${PROJECT_SOURCE_DIR}/src/paramsls.inc
    ${PROJECT_SOURCE_DIR}/src/spin.f)

#add_library(tbpluslibs
#        ${PROJECT_SOURCE_DIR}/src/tbplusexp6all.f
#        ${PROJECT_SOURCE_DIR}/src/tbplusexp6-2.f)

#add_library(mollibs
#        ${PROJECT_SOURCE_DIR}/src/potentials/rohM.f
#        ${PROJECT_SOURCE_DIR}/src/potentials/rohMP.f
#        ${PROJECT_SOURCE_DIR}/src/potentials/rohT.f
#        ${PROJECT_SOURCE_DIR}/src/potentials/roohT.f
#        ${PROJECT_SOURCE_DIR}/src/potentials/ch3ohM.f)

message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

add_executable(auto1dmin.x
    ${PROJECT_SOURCE_DIR}/src/auto1dmin.f)

target_link_libraries(auto1dmin.x 
        ljlibs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS auto1dmin.x DESTINATION bin)
