cmake_minimum_required(VERSION 3.16)

project(VlasovPoissonSolver LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Tell CMake where to find custom Find*.cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- Dependencies -------------------------------------------------------
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(Jansson REQUIRED)  # uses cmake/FindJansson.cmake

# --- Executable ---------------------------------------------------------
add_executable(vlasov_bgk
    src/Vlasov-Poisson_Solver_git.cpp
)

# --- Includes -----------------------------------------------------------
target_include_directories(vlasov_bgk
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${HDF5_INCLUDE_DIRS}
        ${JANSSON_INCLUDE_DIRS}
)

# --- Link libraries -----------------------------------------------------
target_link_libraries(vlasov_bgk
    PRIVATE
        MPI::MPI_CXX
        ${HDF5_LIBRARIES}
        ${JANSSON_LIBRARIES}
)

# Some systems expose an imported hdf5 target; add it if available
if(TARGET hdf5)
    target_link_libraries(vlasov_bgk PRIVATE hdf5)
endif()
