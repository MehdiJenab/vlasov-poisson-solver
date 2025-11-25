# FindJansson.cmake
# Locate Jansson JSON library
#
# This module defines:
#  JANSSON_FOUND - System has Jansson
#  JANSSON_INCLUDE_DIRS - The Jansson include directories
#  JANSSON_LIBRARIES - The libraries needed to use Jansson
#  JANSSON_VERSION - The version of Jansson found

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_JANSSON QUIET jansson)
endif()

# Find the include directory
find_path(JANSSON_INCLUDE_DIR
    NAMES jansson.h
    HINTS
        ${PC_JANSSON_INCLUDEDIR}
        ${PC_JANSSON_INCLUDE_DIRS}
    PATHS
        /usr/include
        /usr/local/include
        /opt/local/include
)

# Find the library
find_library(JANSSON_LIBRARY
    NAMES jansson
    HINTS
        ${PC_JANSSON_LIBDIR}
        ${PC_JANSSON_LIBRARY_DIRS}
    PATHS
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        /usr/lib/x86_64-linux-gnu
)

# Extract version from jansson.h if found
if(JANSSON_INCLUDE_DIR AND EXISTS "${JANSSON_INCLUDE_DIR}/jansson.h")
    file(STRINGS "${JANSSON_INCLUDE_DIR}/jansson.h" JANSSON_VERSION_LINE
         REGEX "^#define[ \t]+JANSSON_VERSION[ \t]+\"[^\"]*\"")
    if(JANSSON_VERSION_LINE)
        string(REGEX REPLACE "^#define[ \t]+JANSSON_VERSION[ \t]+\"([^\"]*)\".*" "\\1"
               JANSSON_VERSION "${JANSSON_VERSION_LINE}")
    endif()
endif()

# Handle standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Jansson
    REQUIRED_VARS
        JANSSON_LIBRARY
        JANSSON_INCLUDE_DIR
    VERSION_VAR
        JANSSON_VERSION
)

# Set output variables
if(JANSSON_FOUND)
    set(JANSSON_LIBRARIES ${JANSSON_LIBRARY})
    set(JANSSON_INCLUDE_DIRS ${JANSSON_INCLUDE_DIR})
    
    # Create imported target
    if(NOT TARGET Jansson::Jansson)
        add_library(Jansson::Jansson UNKNOWN IMPORTED)
        set_target_properties(Jansson::Jansson PROPERTIES
            IMPORTED_LOCATION "${JANSSON_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${JANSSON_INCLUDE_DIR}"
        )
    endif()
endif()

# Mark variables as advanced
mark_as_advanced(
    JANSSON_INCLUDE_DIR
    JANSSON_LIBRARY
)