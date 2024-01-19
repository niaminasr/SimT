set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Source files
set(SOURCES
        PastixSolver.cpp
        PastixSolverParams.cpp
        SolverFabric.cpp
        System.cpp
        )

# Header files
set(HEADERS
        BiConjugateGradientSolver.hpp
        PastixSolver.hpp
        PastixSolverParams.hpp
        Solver.hpp
        SolverFabric.hpp
        SolverParams.hpp
        System.hpp
        )

# Required dependencies
find_package(SPM REQUIRED)
find_package(PASTIX REQUIRED)
find_package(SCOTCH REQUIRED)

# Add the library target
add_library(solver STATIC ${SOURCES} ${HEADERS})

# Set include directories for this component
target_include_directories(solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PASTIX_INC_DIR} ${SCOTCH_INCLUDE_DIR})
target_link_directories(solver PUBLIC ${PASTIX_LIB_DIR} ${SPM_LIB_DIR} ${SCOTCH_LIBRARY_DIR})

# Link against any required libraries
message("MPI:" ${MPI_LIBRARIES})
message("pastix:" ${PASTIX_LIB_DIR} ${PASTIX_LIBRARIES})
target_link_libraries(solver PRIVATE ${MPI_LIBRARIES} spm geometry assemble io tools pastix scotch pastix_kernels spm hwloc lapacke openblas scotcherrexit pthread z m rt ${MPI_LIBRARIES})
