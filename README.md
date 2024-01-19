# EIT Simulation Solver

This is a C++ parallel solver for an EIT (Electrical Impedance Tomography) simulation. It includes modules for assembling, solving, and outputting results.

## Installation

To install, clone the repository and use CMake to build the project. The following dependencies are required:

- [MPI](https://www.open-mpi.org/)
- [PaStiX](https://gitlab.inria.fr/solverstack/pastix)

## Build

To build the project, run the following commands from the root directory of the project:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

The solver is run through the main executable `eit`. 
From the `build` directory, you can run the binary like so:

```bash
    ./eit
```

It is possible to modify the parameters of the simulation without running it again by adding options like so :

```bash
    ./eit -d 100 -c 0.4 -r 1.2   
```
For instance, this runs the simulation with dimensions 100*100, coeff = 0.4, radius = 1.2

To have a complete list of which options are available, run

```bash
    ./eit -h
```

## File Structure

- `src/assemble`: Assembling modules for the simulation
- `src/geometry`: Geometry and grid modules
- `src/io/input`: Configuration file parsing and string tools
- `src/io/output`: Output writers
- `src/solver`: Solver and matrix system modules
- `src/tools`: Miscellaneous tools for the simulation
- `src/main.cpp`: Main executable

## Contributing

Certainly! Here are the instructions for adding a new solver, formatted as a code block:

markdown

### Adding a Solver

To add a new solver to the project, follow these steps:

1. Create a new class that inherits from the `SolverParams` class and name it `<Name>SolverParam`. This class will hold the specific parameters for the new solver. Implement the following functions within this class:

```cpp
class <Name>SolverParam : public SolverParams {
public:
    void setDefault() {
        // Set default values for solver parameters
       // ...
   }

    int getSolverArgc() const {
       // Return the number of arguments required by the solver
        // ...
    }

    char** getSolverArgv() const {
       // Return the array of arguments required by the solver
        // ...
    }
};
```

2. Add a new class attribute `_\<name>SolverParams` of type <Name>SolverParams to the SolverFabric class. This attribute will hold the instance of the solver parameters for the new solver.

```cpp

class SolverFabric {
private:
    // Existing solver attributes
    // ...

    <Name>SolverParams _<name>SolverParams; // Add this line
};
```

Create a new variable <NAME>_SOLVER in the SolverType enum class. This variable will represent the new solver type.

```cpp
enum class SolverType {
    // Existing solver types
    // ...

    <NAME>_SOLVER // Add this line
};
```
Modify the switch cases in the code that handle SolverType variables. Add a new case for the `<NAME>_SOLVER` and handle the specific logic for the new solver.

```cpp
switch (solverType) {
    // Existing solver cases
    // ...

    case SolverType::<NAME>_SOLVER:
        // Handle the logic for the new solver
        // ...
        break;
}
```

Update the -h option in EitOptions to provide information about the value to be used for the new solver. This will help users understand how to select and specify the new solver when running the simulation.

```cpp
        // eitGetOptions switch case
        case 'h':
        default:
            fprintf( stderr,
                     "Usage: -r radius of the form\n"
                     "       -c coeff of the form\n"
                     "       -n number of electrodes\n"
                     "       -d grid dimensions\n"
                     "       -c coefficient\n"
                     "       -s chosen solver   [1 - SolverType::PASTIX_SOLVER]\n"
                     // Other solver types
                     "                          [n - SolverType::<NAME>_SOLVER]\n" );
```

By following these steps, you can successfully add a new solver to the project.

## TO DO

**Add solver:**
- [AMGCL](https://github.com/ddemidov/amgcl)
- [Myramath](https://myramath.org/)