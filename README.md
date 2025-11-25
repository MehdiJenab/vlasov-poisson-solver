# Vlasov–Poisson Kinetic Solver (BGK / Ion-Acoustic / Electron-Hole Dynamics)

This repository contains a modernized CMake-based implementation of a fully kinetic **1D1V Vlasov–Poisson solver** based on the *phase-point (VHS) method*.  
It accurately simulates nonlinear electrostatic structures including:

- BGK modes  
- Ion-acoustic solitons  
- Electron holes and trapped distributions  
- Head-on and overtaking soliton collisions  
- Langmuir-like ionic waves  
- Dust-ion acoustic structures  
- Long-time kinetic evolution with minimal numerical noise  

The code has been used extensively across more than a decade of peer-reviewed publications in **Physics of Plasmas, Physical Review E, IEEE Transactions on Plasma Science**, and **EPJD**.  
A full publication list with DOI links is included below.

---

## 🔥 Features

- Fully kinetic **Vlasov–Poisson** solver (1D1V)
- **VHS / Phase-Point** algorithm (noise-free, high-accuracy)
- **Schamel** trapped electron distribution function (β parameter)
- MPI parallelized Poisson solver
- JSON configuration (via **Jansson**)
- HDF5 output support
- C++17, CMake build system
- Example configurations included

---

## 🛠 Build Instructions

### Requirements
- C++17 compiler  
- CMake ≥ 3.16  
- MPI (OpenMPI / MPICH)  
- HDF5 with MPI support  
- Jansson JSON library  

### Building

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

Executable will be generated as:

```
build/vlasov_bgk
```

---

## ▶ Running a Simulation

From inside `build/`:

```bash
cp ../examples/config*.json .
echo "{}" > Inversed_Matrix.json   # optional placeholder
mpirun -np 4 ./vlasov_bgk
```

---

# 📚 Publications Using This Code

The following peer-reviewed papers were produced directly from this Vlasov–Poisson simulation framework.

---

## **2019 – Physics of Plasmas**

- **Scattering of electron holes in the context of ion-acoustic regime**  
  *Hosseini Jenab, S. M., Spanier, F., Brodin, G.*  
  DOI: https://doi.org/10.1063/1.5055945  

- **Head-on collision of nonlinear solitary solutions to Vlasov–Poisson equations**  
  *Hosseini Jenab, S. M., Brodin, G.*  
  DOI: https://doi.org/10.1063/1.5078865  

---

## **2018 – Physics of Plasmas**

- **A study of the stability properties of Sagdeev solutions in the ion-acoustic regime using kinetic simulations**  
  *Hosseini Jenab, S. M., Spanier, F., Brodin, G.*  
  DOI: https://doi.org/10.1063/1.5036764  

---

## **2017 – IEEE Transactions on Plasma Science / Physical Review E / Physics of Plasmas**

- **Kinetic simulation study of electron holes dynamics during collisions of ion-acoustic solitons**  
  *Hosseini Jenab, S. M., Spanier, F.*  
  DOI: https://doi.org/10.1109/TPS.2017.2715558  

- **Study of ion-acoustic solitons in presence of trapped electrons with a fully kinetic simulation approach**  
  *Hosseini Jenab, S. M., Spanier, F.*  
  DOI: https://doi.org/10.1103/PhysRevE.95.053201  

- **Simulation study of overtaking of ion-acoustic solitons in the fully kinetic regime**  
  *Hosseini Jenab, S. M., Spanier, F.*  
  DOI: https://doi.org/10.1063/1.4978488  

- **Kinetic-simulation study of propagation of Langmuir-like ionic waves in dusty plasma**  
  *Hosseini Jenab, S. M., Spanier, F.*  
  DOI: https://doi.org/10.1109/TPS.2016.2642998  

---

## **2016 – Physics of Plasmas**

- **Study of trapping effect on ion-acoustic solitary waves based on a fully kinetic simulation approach**  
  *Hosseini Jenab, S. M., Spanier, F.*  
  DOI: https://doi.org/10.1063/1.4964909  

---

## **2014 – Physics of Plasmas / European Physical Journal D**

- **Multicomponent kinetic simulation of BGK modes associated with ion acoustic and dust-ion acoustic excitations in electron-ion and dusty plasmas**  
  *Hosseini Jenab, S. M., Kourakis, I.*  
  DOI: https://doi.org/10.1063/1.4869730  

- **Vlasov-kinetic computer simulations of electrostatic waves in dusty plasmas: An overview of recent results**  
  *Hosseini Jenab, S. M., Kourakis, I.*  
  DOI: https://doi.org/10.1140/epjd/e2014-50177-4  

---

## **2011 – Physical Review E / Physics of Plasmas**

- **Preventing the recurrence effect in the Vlasov simulation by randomizing phase-point velocities in phase space**  
  *Abbasi, H., Hosseini Jenab, S. M., Pajouh, H. Hakimi*  
  DOI: https://doi.org/10.1103/PhysRevE.84.036702  

- **Fully kinetic simulation of ion-acoustic and dust-ion acoustic waves**  
  *Hosseini Jenab, S. M., Kourakis, I., Abbasi, H.*  
  DOI: https://doi.org/10.1063/1.3609814  

---

## **2007 – Computer Physics Communications**

- **Vlasov model using kinetic phase point trajectories for the study of BGK modes**  
  *H. Abbasi, M. Ghadimi, Hosseini Jenab, S. M., N. Javaheri*  
  DOI: https://doi.org/10.1016/j.cpc.2007.02.009  

---

# 📂 Repository Structure

```
.
├── cmake/              # CMake module: FindJansson.cmake
├── include/            # Core solver headers
├── src/                # Solver implementation
├── examples/           # Input configuration files
├── CMakeLists.txt      # Build configuration
└── README.md           # Project documentation
```

---

# 🧭 Future Development

- Modern C++ refactoring (RAII, smart pointers, ranges)
- Remove deprecated iterator inheritance
- Modularize the solver into independent components
- Add unit tests (GoogleTest)
- Add GitHub Actions CI
- Python HDF5 visualization tools
- Performance benchmarks (MPI scaling)
