# Source files
set(SOURCES
        EitElectrode.cpp
        EitMath.cpp
        EitVector.cpp
        )

# Header files
set(HEADERS
        EitArray.hpp
        EitBiVector.hpp
        EitElectrode.hpp
        EitMath.hpp
        EitVector.hpp
        )

# Library target
add_library(tools STATIC ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})