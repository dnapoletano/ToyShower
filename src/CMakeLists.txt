add_compile_options(${myCOMPILE_FLAGS})
AUX_SOURCE_DIRECTORY("${PROJECT_SOURCE_DIR}/src" source)
include_directories ("${PROJECT_SOURCE_DIR}/include")

set(RIVET_INCLUDE "${RIVET_PREFIX}/include")
set(RIVET_LIB "${RIVET_PREFIX}/lib")
include_directories (${RIVET_INCLUDE})
include_directories("/usr/local/include")
link_directories(${RIVET_LIB})
link_directories("/usr/local/lib")
link_libraries(Rivet HepMC)


add_executable(${PROJECT_NAME} ${source})
