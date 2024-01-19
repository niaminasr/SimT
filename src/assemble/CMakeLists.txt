# Source files
set(SOURCES
        electrode.cpp
        interface.cpp
        regular.cpp
        write.cpp
        )

# Library target
add_library(assemble STATIC ${SOURCES})

# Include directories
find_package(SPM REQUIRED)
target_link_directories(assemble PUBLIC ${SPM_LIB_DIR})
target_link_libraries(assemble geometry spm)
target_include_directories(assemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/.. ${SPM_INC_DIR})
