# Source files
set(SOURCES
        input/parameters.cpp
        input/stringTools.cpp
        output/Writers.cpp
        )

# Header files
set(HEADERS
        input/parameters.hpp
        input/stringTools.hpp
        output/Writers.hpp
        )

# Library target
add_library(io STATIC ${SOURCES} ${HEADERS})

# Include directories
target_link_libraries(io geometry)
target_include_directories(io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
