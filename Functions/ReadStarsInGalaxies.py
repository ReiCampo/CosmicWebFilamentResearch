###########################################################################
###########################################################################
###                                                                     ###
###                  PULLING OUT VALUES NEEDED FOR SFH                  ###
###                                                                     ###
###########################################################################
###########################################################################

### Fixes for the code:
# First, use the center of mass of the galaxy to recenter the the positions on
# the galaxy

# You also do it for the velocities too. 

# Put yourself in the frame of reference of the galaxy, once you do that then 
# you can start to calculate velocity dispersions. This will make your graph
# make more sense!

# When there is a huge gap in your positions, there is perodicity that you need
# to consider. 

# We only identify the zoom regions in New Horizons.

# You can look at v_rotational on the line of sight
# YOu can choose any axis to be your line of sight. 

# You can read the angular momentum of the galaxy to determine how the galaxy
# is tilted. You can project on your new x prime and y prime axis to see how the
# galaxy is oriented.


# Add into compute stars to do edge on/recentering of the stars


##----------------------------------------------------------------
##                  Importing Necessary Packages                 -
##----------------------------------------------------------------

import os

# This sets the working directory to the folder where this script lives. Make 
# sure that your "read_treebricks" file is in the same location as this script!
os.chdir("/Users/RachelCampo/Desktop/CUNY/Filament Research/CosmicWebFilamentResearch/Functions")

# This is for active development of the treebricks_function and compute_stars 
# code:
import importlib
import treebricks_function
import compute_stars
importlib.reload(treebricks_function)
importlib.reload(compute_stars)

# Adding in Treebrick functions:
from treebricks_function import ReadHaloShapes, ReadTreebrick_highp, ReadTreebrick_lowp, GalaxyCatalog
from compute_stars import ReadGalData, RecenterStars

from os.path import join
import re
import numpy as np
import time
import struct
import matplotlib.pyplot as plt
import pandas as pd

# Packages to help with translation from Fortran to Python:
import subprocess
from scipy.io import FortranFile
import ctypes as c


##----------------------------------------------------------------
##                        Importing Data                         -
##----------------------------------------------------------------


data_path = "/Users/RachelCampo/Desktop/CUNY/Filament Research/Star Formation Rate In Filaments/New Horizons Data/TreeStarsAdaptaHOP/"

# This data came from:
# /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/GAL_00970 on Infinity
galaxy970_data = data_path + str("GAL_00970/gal_stars_0004071")

# This data came from:
# /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/GAL_00100 on Infinity
galaxy100_data = data_path + str("GAL_00100/gal_stars_0001297")


##----------------------------------------------------------------
##                          Reading Data                         -
##----------------------------------------------------------------

obj970 = ReadGalData(gal_data = galaxy970_data)
stars_gals970_df = obj970.star_info
# Now I am going to use the new function, RecenterStars:

recentered_stars_gals970_df = RecenterStars(gal_info = obj970.galaxy_header_info,
                                            star_info = obj970.star_info)
adjusted_stars_970_df = recentered_stars_gals970_df.adjusted_star_info
print(adjusted_stars_970_df)

# obj100 = ReadGalData(gal_data = galaxy100_data)
# stars_gals100_df = obj100.galaxy_info
# print(stars_gals100_df)


###  Now I am going to calculate some values for plotting later. I will start   
###  by calculating the total radius and velocity. 

stars_gals970_df = adjusted_stars_970_df.dropna(subset=["Star_Position_X (Mpc)", 
                                                        "Star_Position_Y (Mpc)",
                                                        "Star_Position_Z (Mpc)",
                                                        "Star_Velocity_X (km/s)",
                                                        "Star_Velocity_Y (km/s)",
                                                        "Star_Velocity_Z (km/s)"])

# First, I'm pulling out the positions and velocities to use later:
x_pos_970 = stars_gals970_df["Star_Position_X (Mpc)"]
y_pos_970 = stars_gals970_df["Star_Position_Y (Mpc)"]
z_pos_970 = stars_gals970_df["Star_Position_Z (Mpc)"]   

x_vel_970 = stars_gals970_df["Star_Velocity_X (km/s)"]
y_vel_970 = stars_gals970_df["Star_Velocity_Y (km/s)"]
z_vel_970 = stars_gals970_df["Star_Velocity_Z (km/s)"]                          

# Doing some inspecting, I've noticed that there are some stars that are way 
# beyond a physical galaxies distance.

# Now I am calculating the total radius and velocities:
stars_970_total_radius = np.sqrt(x_pos_970 ** 2 + y_pos_970 ** 2 + z_pos_970 ** 2)

stars_970_total_velocity = np.sqrt(x_vel_970 ** 2 + y_vel_970 ** 2 + z_vel_970 ** 2)

#star_mask = (np.isfinite(stars_970_total_velocity) & (stars_970_total_radius >= 3.5) & (stars_970_total_radius < 3.8))

# Apply the mask:
#stars_970_total_radius_adjusted = stars_970_total_radius[star_mask]

#stars_970_total_velocity_adjusted = stars_970_total_velocity[star_mask]

# I am finding the minimum and maximum radius to use for binning the radius
# later:
min_rad_970 = stars_970_total_radius.min()

max_rad_970 = stars_970_total_radius.max()

bin_number = 100

# Setting up bins for the plot. I'm going to start with 50 bins for now:
radius_bins_970 = np.linspace(min_rad_970, max_rad_970, bin_number + 1)

bin_width = radius_bins_970[1] - radius_bins_970[0]

print(min_rad_970)
print(max_rad_970)
print(bin_width)

# Setting up lists to use later for plotting:
centered_radius = []
velocity_dispersion = []
dispersion_error = []


for i in radius_bins_970:
    
    bin_mask = (stars_970_total_radius_adjusted >= i) & (stars_970_total_radius_adjusted < i + bin_width)
    
    stars_in_bin = stars_970_total_velocity_adjusted[bin_mask]
    
    print(f"Bin {i:.3f} to {i + bin_width:.3f}: {len(stars_in_bin)} stars")
    
    bin_centers = i + bin_width / 2
    centered_radius.append(bin_centers)
    
    sigma = stars_in_bin.std(ddof = 1)
    velocity_dispersion.append(sigma)
    
    sigma_error = sigma / np.sqrt(2 * (len(stars_in_bin) - 1))
    dispersion_error.append(sigma_error)
    

##----------------------------------------------------------------
##                            Plotting                           -
##----------------------------------------------------------------

fig, ax = plt.subplots()

ax.errorbar(x = centered_radius,
            y = velocity_dispersion,
            yerr = dispersion_error,
            fmt = "o",
            color = "#2ABF83",
            ecolor = "gray",
            markeredgecolor = "black",
            markersize = 10,
            label = f"Velocity Dispersion, Bin = {bin_number}")
ax.set_xscale("log")
ax.set_title("The Velocity Dispersion for Galaxy 0004071 at Time Step 970")
ax.set_xlabel("Bin Radii in Log Scale (Mpc)")
ax.set_ylabel("Velocity Dispersion at Bin with Errors (km/s)")
ax.legend()
ax.grid(True, alpha = 0.3)

plt.show()


plt.plot(np.log10(x_pos_970), np.log10(y_pos_970), color = "green")
plt.show()