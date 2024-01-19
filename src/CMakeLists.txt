set(THIS_TARGET eit)

# Add source files
set(SOURCES
        main.cpp
        EitOptions.cpp
        )

# Add subdirectories for different components
add_subdirectory(assemble)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(solver)
add_subdirectory(tools)

# Add the main executable
add_executable(${THIS_TARGET} ${SOURCES})

# Link against the required libraries
target_link_libraries(${THIS_TARGET} PUBLIC assemble geometry io solver tools ${MPI_LIBRARIES})
set_target_properties(${THIS_TARGET} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(${THIS_TARGET} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")