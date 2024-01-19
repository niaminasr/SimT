# Source files
set(SOURCES
        Form.cpp
        Grid.cpp
        Point.cpp
        )

# Header files
set(HEADERS
        Form.hpp
        Grid.hpp
        Point.hpp
        )

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Library target
add_library(geometry STATIC ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link
target_link_libraries(geometry PUBLIC tools ${MPI_LIBRARIES})
