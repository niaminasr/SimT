set(THIS_TARGET tests)


find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(tests_point tests_point.cpp)
add_executable(tests_form tests_form.cpp)
add_executable(tests_assemble tests_assemble.cpp)

target_link_libraries(tests_point PRIVATE ${MPI_LIBRARIES} ${EIGEN3_LIBRARIES} assemble geometry io solver tools)
set_target_properties(tests_point PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(tests_point PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}") 

target_link_libraries(tests_form ${MPI_LIBRARIES} ${EIGEN3_LIBRARIES} assemble geometry io solver tools)
set_target_properties(tests_form PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(tests_form PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}") 

target_link_libraries(tests_assemble PRIVATE ${MPI_LIBRARIES} ${EIGEN3_LIBRARIES} assemble geometry io solver tools)
set_target_properties(tests_assemble PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(tests_assemble PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

##Test Point.cpp
add_test(NAME Point_getElectrodeIdForAngle COMMAND tests_point 1)
add_test(NAME Point_computeGlobalIndex COMMAND tests_point 2)
add_test(NAME Point_setAsRegular COMMAND tests_point 3)
add_test(NAME Point_setAsInterface COMMAND tests_point 4)

##Test Form.cpp
add_test(NAME Form_radius COMMAND tests_form 1)
add_test(NAME Form_firstDerivativeRadius COMMAND tests_form 2)
add_test(NAME Form_secondDerivativeRadius COMMAND tests_form 3)
add_test(NAME Form_getEndAngleForLength COMMAND tests_form 4)
add_test(NAME Form_buildThetas COMMAND tests_form 5)
add_test(NAME Form_derivativeLevelSet COMMAND tests_form 6)
add_test(NAME Form_computeInterfacePoint COMMAND tests_form 7)

#Test Electrode.cpp
add_test(NAME Full_assemble COMMAND tests_assemble 1)
