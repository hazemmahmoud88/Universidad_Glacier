"""
Nuth and Kääb coregistration
============================

Nuth and Kääb (`2011 <https://doi.org/10.5194/tc-5-271-2011>`_) coregistration allows for horizontal and vertical shifts to be estimated and corrected for.
In this script, this approach is implemented through a coregistration class.

For more information about the approach, see the relevant documentation.
"""
import geoutils as gu
import numpy as np

import xdem

# Example files
reference_dem = xdem.DEM(r"path/to/reference_dem.tif")
dem_to_be_aligned = xdem.DEM(r"path/to/dem_to_be_aligned.tif")
mask_out_stable = gu.Vector(r"path/to/stable_mask.shp")

# Create a stable ground mask (not glacierized) to mark "inlier data"
inlier_mask = ~mask_out_stable.create_mask(reference_dem)

# Visualize initial elevation change between reference and aligned DEMs
diff_before = reference_dem - dem_to_be_aligned
diff_before.show(cmap="coolwarm_r", vmin=-100, vmax=100, cb_title="Elevation change (m)")
rmse_b01 = (reference_dem - dem_to_be_aligned)**2

# Estimate and apply horizontal and vertical shifts using coregistration
coreg = xdem.coreg.NuthKaab()
coreg.fit(reference_dem.data, dem_to_be_aligned.data, inlier_mask, transform=reference_dem.transform)
aligned_dem_data = coreg.apply(dem_to_be_aligned)
aligned_dem_data.save(r"path/to/save/aligned_dem.tif")

# Visualize improved elevation change after coregistration
diff_after = reference_dem - aligned_dem_data
diff_after.show(cmap="coolwarm_r", vmin=-100, vmax=100, cb_title="Elevation change (m)")
rmse_a01 = (reference_dem - dem_to_be_aligned)**2

# Evaluate error metrics to validate improvement
inliers_before = diff_before.data[inlier_mask].compressed()
med_before, nmad_before = np.median(inliers_before), xdem.spatialstats.nmad(inliers_before)
diff_mean_b = np.mean(inliers_before)
rmse_b02 = rmse_b01.data[inlier_mask].compressed()
rmse_b03 = np.sqrt(np.mean(rmse_b02))

inliers_after = diff_after.data[inlier_mask].compressed()
med_after, nmad_after = np.median(inliers_after), xdem.spatialstats.nmad(inliers_after)
diff_mean_a = np.mean(inliers_after)
rmse_a02 = rmse_a01.data[inlier_mask].compressed()
rmse_a03 = np.sqrt(np.mean(rmse_a02))

print(f"Error before coregistration: median = {med_before:.2f} - NMAD = {nmad_before:.2f} m")
print(f"Error after coregistration: median = {med_after:.2f} - NMAD = {nmad_after:.2f} m")
print(f"Mean error before coregistration: mean = {diff_mean_b:.2f} m")
print(f"Mean error after coregistration: mean = {diff_mean_a:.2f} m")
print(f"RMSE before coregistration: RMSE = {rmse_b03:.2f} m")
print(f"RMSE after coregistration: RMSE = {rmse_a03:.2f} m")

# Iteratively apply the coregistration to multiple DEMs
for i in range(1, 10):
    aligned_file = f"aligned_dem_{i}.tif"
    print(aligned_file)
    dem_to_be_aligned = xdem.DEM(f"path/to/dem_folder/{aligned_file}")

    diff_before = reference_dem - dem_to_be_aligned
    rmse_b01 = (reference_dem - dem_to_be_aligned) ** 2

    coreg = xdem.coreg.NuthKaab()
    coreg.fit(reference_dem.data, dem_to_be_aligned.data, inlier_mask, transform=reference_dem.transform)
    aligned_dem_data = coreg.apply(dem_to_be_aligned)
    next_aligned_file = f"aligned_dem_{i + 1}.tif"
    print(next_aligned_file)
    aligned_dem_data.save(f"path/to/save/{next_aligned_file}")

    diff_after = reference_dem - aligned_dem_data
    rmse_a01 = (reference_dem - dem_to_be_aligned) ** 2

    inliers_before = diff_before.data[inlier_mask].compressed()
    med_before, nmad_before = np.median(inliers_before), xdem.spatialstats.nmad(inliers_before)
    diff_mean_b = np.mean(inliers_before)
    rmse_b02 = rmse_b01.data[inlier_mask].compressed()
    rmse_b03 = np.sqrt(np.mean(rmse_b02))

    inliers_after = diff_after.data[inlier_mask].compressed()
    med_after, nmad_after = np.median(inliers_after), xdem.spatialstats.nmad(inliers_after)
    diff_mean_a = np.mean(inliers_after)
    rmse_a02 = rmse_a01.data[inlier_mask].compressed()
    rmse_a03 = np.sqrt(np.mean(rmse_a02))

    print(f"Error before coregistration: median = {med_before:.2f} - NMAD = {nmad_before:.2f} m")
    print(f"Error after coregistration: median = {med_after:.2f} - NMAD = {nmad_after:.2f} m")
    print(f"Mean error before coregistration: mean = {diff_mean_b:.2f} m")
    print(f"Mean error after coregistration: mean = {diff_mean_a:.2f} m")
    print(f"RMSE before coregistration: RMSE = {rmse_b03:.2f} m")
    print(f"RMSE after coregistration: RMSE = {rmse_a03:.2f} m")

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

# Load the final aligned DEM and the original reference DEM
aligned_dem = gdal.Open(r"path/to/final_aligned_dem.tif")
reference_clip = gdal.Open(r"path/to/reference_dem.tif")

# Clip the reference and aligned DEMs to the same extent using a shapefile mask
ref_clip = gdal.Warp("path/to/pro_ref_dem_clip.tif", reference_clip,
                     cutlineDSName="path/to/clip_shapefile.shp", cropToCutline=True, dstNodata=np.nan)
ali_clip = gdal.Warp("path/to/pro_aligned_dem_clip.tif", aligned_dem,
                     cutlineDSName="path/to/clip_shapefile.shp", cropToCutline=True, dstNodata=np.nan)

# Load the clipped DEMs into xdem
ref_clip_dem = xdem.DEM("path/to/pro_ref_dem_clip.tif")
ali_clip_dem = xdem.DEM("path/to/pro_aligned_dem_clip.tif")

# Calculate and plot the difference between the clipped DEMs
diff_after = ref_clip_dem - ali_clip_dem
diff_after.show(cmap="coolwarm_r", vmin=-100, vmax=100, cb_title="Elevation change (m)")
rmse_b01 = (ref_clip_dem - ali_clip_dem) ** 2

# Calculate and print final error statistics
inliers_after = diff_after.data
med_after, nmad_after = np.median(inliers_after), xdem.spatialstats.nmad(inliers_after)
diff_mean_a = np.mean(inliers_after)
rmse_a02 = rmse_b01.data
rmse_a03 = np.sqrt(np.mean(rmse_a02))

print(f"Final Error: median = {med_after:.2f} - NMAD = {nmad_after:.2f} m")
print(f"Final Mean error: mean = {diff_mean_a:.2f} m")
print(f"Final RMSE: RMSE = {rmse_a03:.2f} m")
