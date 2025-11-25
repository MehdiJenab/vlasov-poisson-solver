# Vlasov–Poisson Kinetic Solver

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/CMake-3.16+-blue.svg)](https://cmake.org/)

A fully kinetic **1D1V Vlasov–Poisson solver** implementing the phase-point (VHS) method for simulating nonlinear electrostatic structures in plasmas. This code has been used in **15+ peer-reviewed publications** spanning over a decade of research in plasma physics.

## Overview

This solver accurately captures kinetic physics of:

- **BGK modes** and electron phase-space holes
- **Ion-acoustic solitons** with trapped particle distributions
- **Soliton collisions** (head-on and overtaking)
- **Langmuir-like ionic waves** in dusty plasmas
- **Dust-ion acoustic structures**
- Long-time evolution with minimal numerical noise

The phase-point method tracks individual phase-space trajectories, enabling noise-free simulations of kinetic plasma phenomena that are difficult to capture with traditional PIC codes.

## Key Features

- **Fully kinetic**: No approximations to the Vlasov equation
- **Phase-point (VHS) algorithm**: High accuracy with minimal numerical diffusion
- **Schamel distribution**: Built-in support for trapped electron populations (β parameter)
- **MPI parallelization**: Scalable across multiple nodes
- **JSON configuration**: Easy parameter management via Jansson library
- **HDF5 output**: Standard format for analysis and visualization
- **Modern C++17**: Clean, maintainable codebase
- **CMake build system**: Cross-platform compatibility

## Quick Start

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+)
- CMake ≥ 3.16
- MPI implementation (OpenMPI, MPICH, Intel MPI)
- HDF5 with parallel I/O support
- Jansson JSON library

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install cmake build-essential
sudo apt-get install libopenmpi-dev libhdf5-openmpi-dev libjansson-dev
```

#### macOS (Homebrew)
```bash
brew install cmake open-mpi hdf5-mpi jansson
```

#### CentOS/RHEL
```bash
sudo yum install cmake gcc-c++
sudo yum install openmpi-devel hdf5-openmpi-devel jansson-devel
```

### Building

```bash
git clone https://github.com/MehdiJenab/vlasov-poisson-solver.git
cd vlasov-poisson-solver
mkdir build && cd build
cmake ..
cmake --build . -j4
```

The executable `vlasov_bgk` will be created in the `build/` directory.

### Running Your First Simulation

```bash
# From the build directory
cp ../examples/config.json .
mpirun -np 4 ./vlasov_bgk config.json
```

Output files in HDF5 format will be written to the current directory.

## Example Configurations

Three example configurations are provided in `examples/`:

- **`config.json`**: Basic BGK mode simulation
- **`config_solitons.json`**: Ion-acoustic soliton propagation
- **`config_output.json`**: Detailed output configuration example

### Configuration Parameters

Key parameters you can adjust:

```json
{
  "grid": {
    "nx": 1024,           // Spatial grid points
    "nv": 256,            // Velocity grid points
    "xmin": 0.0,
    "xmax": 100.0,
    "vmin": -6.0,
    "vmax": 6.0
  },
  "time": {
    "dt": 0.01,           // Time step
    "nsteps": 10000       // Total steps
  },
  "species": [
    {
      "name": "electrons",
      "mass": 1.0,
      "charge": -1.0,
      "temperature": 1.0,
      "beta": 0.2         // Trapped fraction (Schamel)
    }
  ]
}
```

See the example files for complete configuration options.

## Physical Models

### Schamel Distribution

The solver supports the Schamel distribution function for trapped particles:

```
f(v) = A · exp(-ε/T) / (1 + β·sinh²(√(ε/T)))
```

where β controls the trapped particle population, enabling accurate modeling of BGK modes and electron holes.

### Poisson Solver

The electrostatic field is computed via MPI-parallelized solution of:

```
∇²φ = -ρ/ε₀
E = -∇φ
```

## Output Format

Simulation data is written in HDF5 format with the following structure:

```
output.h5
├── timestep_0000/
│   ├── distribution_electrons
│   ├── distribution_ions
│   ├── electric_field
│   └── potential
├── timestep_0001/
│   └── ...
└── metadata
    ├── grid_parameters
    └── physical_constants
```

## Publications

This code has been used in the following peer-reviewed research:

### Recent (2017-2019)

1. **Hosseini Jenab, S. M., Spanier, F., Brodin, G.** (2019)  
   *Scattering of electron holes in the context of ion-acoustic regime*  
   Physics of Plasmas **26**, 022305 | [DOI: 10.1063/1.5055945](https://doi.org/10.1063/1.5055945)

2. **Hosseini Jenab, S. M., Brodin, G.** (2019)  
   *Head-on collision of nonlinear solitary solutions to Vlasov–Poisson equations*  
   Physics of Plasmas **26**, 012107 | [DOI: 10.1063/1.5078865](https://doi.org/10.1063/1.5078865)

3. **Hosseini Jenab, S. M., Spanier, F., Brodin, G.** (2018)  
   *Stability properties of Sagdeev solutions in the ion-acoustic regime*  
   Physics of Plasmas **25**, 062307 | [DOI: 10.1063/1.5036764](https://doi.org/10.1063/1.5036764)

4. **Hosseini Jenab, S. M., Spanier, F.** (2017)  
   *Electron holes dynamics during collisions of ion-acoustic solitons*  
   IEEE Trans. Plasma Sci. **45**(8), 2022 | [DOI: 10.1109/TPS.2017.2715558](https://doi.org/10.1109/TPS.2017.2715558)

5. **Hosseini Jenab, S. M., Spanier, F.** (2017)  
   *Ion-acoustic solitons with trapped electrons: Fully kinetic approach*  
   Phys. Rev. E **95**, 053201 | [DOI: 10.1103/PhysRevE.95.053201](https://doi.org/10.1103/PhysRevE.95.053201)

### Earlier Work (2007-2016)

<details>
<summary>Click to expand full publication list</summary>

6. **Hosseini Jenab, S. M., Spanier, F.** (2017)  
   *Overtaking of ion-acoustic solitons in the fully kinetic regime*  
   Physics of Plasmas **24**, 032309 | [DOI: 10.1063/1.4978488](https://doi.org/10.1063/1.4978488)

7. **Hosseini Jenab, S. M., Spanier, F.** (2017)  
   *Langmuir-like ionic waves in dusty plasma: Kinetic simulation*  
   IEEE Trans. Plasma Sci. **45**(3), 413 | [DOI: 10.1109/TPS.2016.2642998](https://doi.org/10.1109/TPS.2016.2642998)

8. **Hosseini Jenab, S. M., Spanier, F.** (2016)  
   *Trapping effect on ion-acoustic solitary waves: Fully kinetic approach*  
   Physics of Plasmas **23**, 082104 | [DOI: 10.1063/1.4964909](https://doi.org/10.1063/1.4964909)

9. **Hosseini Jenab, S. M., Kourakis, I.** (2014)  
   *Multicomponent kinetic simulation of BGK modes in dusty plasmas*  
   Physics of Plasmas **21**, 032112 | [DOI: 10.1063/1.4869730](https://doi.org/10.1063/1.4869730)

10. **Hosseini Jenab, S. M., Kourakis, I.** (2014)  
    *Vlasov-kinetic simulations of electrostatic waves in dusty plasmas*  
    Eur. Phys. J. D **68**, 77 | [DOI: 10.1140/epjd/e2014-50177-4](https://doi.org/10.1140/epjd/e2014-50177-4)

11. **Abbasi, H., Hosseini Jenab, S. M., Pajouh, H. Hakimi** (2011)  
    *Preventing recurrence effect in Vlasov simulation*  
    Phys. Rev. E **84**, 036702 | [DOI: 10.1103/PhysRevE.84.036702](https://doi.org/10.1103/PhysRevE.84.036702)

12. **Hosseini Jenab, S. M., Kourakis, I., Abbasi, H.** (2011)  
    *Fully kinetic simulation of ion-acoustic and dust-ion acoustic w