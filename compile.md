## Required Modules

```sh
module load cmake
module load gcc/11.3.0
module unload openmpi
module load intel-mpi/2021.9.0
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc12-impi
export NETCDF_LIBRARIES=$(nc-config --libs)
export NETCDF_INCLUDES=$(nc-config --includedir)
export CFLAGS="-std=c11"
export CXXFLAGS="-std=c++20"
```

## To Compile

```sh
cmake ..
make
```

## If You Have Compilation Problems

```sh
make clean_all
```