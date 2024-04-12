# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
ds_file =  snakemake.input[0]
area_file =  snakemake.input[2]
out_file = snakemake.output[0]
      
# test
ds_file = "results/downscaled/Côte-d'Ivoire_CMIP6_world_MIROC_MIROC-ES2L_ssp126_r8i1p1f2_none_none_chelsa2_monthly-means_2006-2019_1980-2005_bc.nc"
area_file = "results/areas/Côte-d'Ivoire.shp"

# libs
import xarray as xr
import geopandas as gp
import numpy as np

# code
tmin = 24
prmin = 2000
area_shp = gp.read_file(area_file)
ds = xr.open_dataset(ds_file).rio.clip(area_shp.geometry.values, area_shp.crs)
ds_year = ds[["pr", "tas"]].mean("month")
ds_year['pr'] = ds[["pr"]].sum("month").pr
ds_year = ds_year.assign(tmf = lambda x: np.logical_and((x.tas > 20), (x.pr > 1500)))


# view
import matplotlib.pyplot as plt
ds_year.tmf.plot()
plt.show()
