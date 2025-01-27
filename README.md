SWE
===


## About the project

SWE are shallow water equation solvers, and we have focused on optimizing the code for performance. We did resolving of load imbalances, vectorization and parallelization.


## Compilation
Firstly we need to get the git submodles: 

```sh
git submodule init
git submodule update
```

To compile the code you need the following modules:

```sh
module load cmake
module load gcc/11.3.0
module unload openmpi
module load intel-mpi/2021.9.0
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc12-impi
```

In addition, you need to set the following environment variables:

```sh
export NETCDF_LIBRARIES=$(nc-config --libs)
export NETCDF_INCLUDES=$(nc-config --includedir)
export CFLAGS="-std=c11"
export CXXFLAGS="-std=c++20"
```

Then create a build directory and navigate into it:

```sh
mkdir build
cd build
```

Lastly, run CMake and Make:

```sh
cmake ..
make
```

If you run into errors, try to delete the content of the build folder, and run cmake again:

## Running the code

* Run the code in serial via `./SWE-MPI-Runner`
* Run the code in parallel via `mpirun -np nproc ./SWE-MPI-Runner`
* With `./SWE-MPI-Runner --help`, you can see additional command-line arguments you can pass.

For example, a command to run the code in parallel with 10 processes would be:

```sh
mpiexec -np 10 ./SWE-MPI-Runner -x 500 -y 500
```

## Visualize the Results

The command line version of SWE will write a netCDF file or multiple ASCII-VTK files (depending on the build configuration) which can be opened and visualized with ParaView.

**Hint:** The output files contain the absolute height of the water column _h_. To get the relative height of the water column, use a _Calculator_ in ParaView with the formula `h+b`. If you have dry cells in your scenario, you may want to use the formula `min(h, h+b)`. This will give you the relative height for wet cells and 0 for dry cells.