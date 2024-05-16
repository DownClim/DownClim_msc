# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
ds_file =  snakemake.input[0]
rsds_file = snakemake.input[1]
area_file =  snakemake.input[2]
out_file = snakemake.output[0]
      
# test
# rsds_file = "results/tmf/rsds/Côte-d'Ivoire.nc"
ds_file = "results/downscaled/Côte-d'Ivoire_CMIP6_world_MIROC_MIROC-ES2L_ssp126_r8i1p1f2_none_none_chelsa2_monthly-means_2006-2019_1980-2005_bc.nc"
area_file = "results/areas/Côte-d'Ivoire.shp"

# libs
import xarray as xr
import geopandas as gp
import numpy as np

# limits
pet_pr_max = 1
tmin = 16
prmin = 1000

# code
area_shp = gp.read_file(area_file)
ds = xr.open_dataset(ds_file).rio.clip(area_shp.geometry.values, area_shp.crs)
ds_rsds = xr.open_dataset(rsds_file).rio.clip(area_shp.geometry.values, area_shp.crs)
ds = xr.merge([ds, ds_rsds])

tmean = ds["tas"] - 273
lat = ds.lat * np.pi/180
pet_oudin = pyet.oudin(tmean, lat=lat)

ds['pet'] = 0.0023 * ds.rsds * (ds.tas + 17.8) * pow((ds.tasmax - ds.tasmin), 0.5)
# using pyet: https://pyet.readthedocs.io/en/dev/examples/09_CMIP6_data.html 
ds.pet.attrs = {'standard_name': 'potential evapotranspiration', 
                'long_name': 'Monthly potential evapotranspiration',
                'units': 'mm month-1', 
                'explanation' : 'Potential evapotranspiration for each month; calculated with the Penman-Monteith equation.'}
ds_year = ds.mean("month")
ds_year['pr'] = ds[["pr"]].sum("month").pr
ds_year['pet'] = ds[["pet"]].sum("month").pet
ds_year['pet_pr'] = ds_year.pet / ds_year.pr
ds_year['tmf'] = (ds_year.tas >= tmin).astype(int)*(ds_year.pr >= prmin).astype(int)*(ds_year.pet_pr <= pet_pr_max).astype(int)
ds_year[['tmf']].to_netcdf(out_file)

# view
# import matplotlib.pyplot as plt
# ds_year.tmf.plot()
# plt.show()