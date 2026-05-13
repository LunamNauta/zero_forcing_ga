# Zero Forcing Genetic Algorithm Solver

This project implements a genetic algorithm to find the minimum Zero Forcing set for a given graph. The end goal is to combine a genetic algorithm based heuristic with an exact solver such as the Fort Cover Model IP to solve for the Zero Forcing number of a graph.

## Project Structure

* include/: Header files
* src/: Source files
* examples/: Example scripts verifying the validity of the algorithms
* build/: Output location for built project

## Build

### Note
The current ```CMakeLists.txt``` assumes that
```
GUROBI_DIR=../gurobi/gurobi1103/linux64
GUROBI_INC=$GUROBI_DIR/include
GUROBI_LIB=$GUROBI_DIR/lib
```
relative to the project root. If this is not the case for you, ensure that these variabes are defined when you call cmake.

### Default
```bash
mkdir -p build && cd build \
&& cmake ../               \
&& make
```

### Minimum Overload (if needed)
```bash
mkdir -p build && cd build                 \
&& cmake ../ -DGUROBI_DIR=<PATH_TO_GUROBI> \
&& make
```