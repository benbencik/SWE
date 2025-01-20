import netCDF4 as nc
import numpy as np

def read_netcdf(file_path):
    dataset = nc.Dataset(file_path, 'r')
    h = dataset.variables['h'][:]
    hu = dataset.variables['hu'][:]
    hv = dataset.variables['hv'][:]
    dataset.close()
    return h, hu, hv

def compare_results(h1, hu1, hv1, h2, hu2, hv2):
    assert np.allclose(h1, h2), "Mismatch in h values"
    assert np.allclose(hu1, hu2), "Mismatch in hu values"
    assert np.allclose(hv1, hv2), "Mismatch in hv values"
    print("All values match!")

if __name__ == "__main__":
    baseline = 'SWE_vecsolver_250x250.nc'
    
    #! change the following line
    test = 'SWE_00.nc' 

    h1, hu1, hv1 = read_netcdf(baseline)
    h2, hu2, hv2 = read_netcdf(test)

    compare_results(h1, hu1, hv1, h2, hu2, hv2)
    