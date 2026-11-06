# This script was written by Rei Campo
# Last updated on 9/15/2026

############################################################################
############################################################################
###                                                                      ###
###                          IMPORTING PACKAGES                          ###
###                                                                      ###
############################################################################
############################################################################



import numpy as np
import pandas as pd
import glob
import argparse
import re
import math
from collections import defaultdict
from scipy.spatial import KDTree

# These packages were created by the CosmicWeb group and serve as helper
# functions:

import filament_tools as ft
from vector import segment_distance


# These packages below were created by Rei Campo to assist in reading
# galaxy treebricks and calculating certain metrics:

from compute_stars import ReadGalData
from calculate_metrics import CalculateMetrics




###########################################################################
###########################################################################
###                                                                     ###
###                   CREATING COMMAND LINE ARGUMENTS                   ###
###                                                                     ###
###########################################################################
###########################################################################



class UserInputs:
    
    def __init__(self):
        
        '''
        This is the initializer function that takes in all user-inputted 
        settings.
            
        ----------------------------------------------------------------------
            
        Inputs:
            
            self.parser (parser):
                The parser that gets created to store all inputs.
                
            self.nh_time_step_inputs (class function):
                This is the function that takes in NewHorizon specific inputs.
                
            self.split_data_inputs (class function):
                This is the function that takes in all other inputs like how
                the data should be split up, how many k-nearest-neighbors and 
                interpolation points to use. 
                
            self.args (parser):
                This is the variable that is going to store all the user-
                inputted values.
                
            self.print_user_outputs (class function):
                This function serves as a readability function to print out 
                all the user-inputted values.
                
        '''
        
        self.parser = argparse.ArgumentParser(description = "Calculating SFH and SFR for Galaxies in the NewHorizon simulation")
        
        self.nh_time_step_inputs()
        self.split_data_inputs()
        
        self.args = self.parser.parse_args()
        
        self.print_user_outputs()



    ##---------------------------------------------------------------
    ##          Inputs on the Time Step and its Cosmology           -
    ##---------------------------------------------------------------


    
    def nh_time_step_inputs(self):
        '''
        This function takes in user inputted data based on the time step the user
        wants to look at and that time step's cosmology. '
        
        The cosmological parameters are all defaulted to time step 1026 in 
        NewHorizon.
        '''
        
        parser = self.parser
        
        
        time_step_flags = {
            
            # The data path for where your galaxies are located:
            "--data_path": {"type": str, "required": True},
            
            # The data path where your filament file is located:
            "--filament_file_path": {"type": str, "required": True},
            
            # The output path where you want the saved data to be housed:
            "--output_path": {"type": str, "required": True},
            
            # This is the selected timestep one would like to explore. The default
            # is set to the latest time step, 1026. All cosmological parameters
            # are defaulted to time step 1026:
            "--time_step": {"type": int, "default": 1026},
            
            # This bin width will be used to determine the SFH and SFR. The default
            # is set to 100 Myrs
            "--bin_width": {"type": int, "default": 100e6},
            
            
            # Below are the arguments needed for the cosmology of the specific
            # time step. You can get this information from the time step's
            # info.txt file:
            "--alpha":     {"type": float, "default": 1e-6},
            "--omega_m":   {"type": float, "default": 0.272000014781952e0},
            "--omega_l":   {"type": float, "default": 0.727999985218048e0},
            "--omega_k":   {"type": float, "default": 0.0e0},
            "--aexp":      {"type": float, "default": 0.853613930315508e0},
            "--ntable":    {"type": int,   "default": 1000},
            "--H0":        {"type": float, "default": 0.704000015258789e2},
            "--omega_r":   {"type": float, "default": 0.0}, # Still have to look into this, but omega_r is negligible
            "--unit_l":    {"type": float, "default": 0.373456098160067E+27},
            
            # If you created a gas cube, these are the values found in your
            # submit script:
            "--boxlen":    {"type": float, "default": 0.100000000000000E+01},
            "--xmax":      {"type": float, "required": True},
            "--xmin":      {"type": float, "required": True}
            }
            
            
            
        # Adding in NH inputs into the parser:    
        for flag, kwargs in time_step_flags.items():
            parser.add_argument(flag, **kwargs)
    
    
    
    
    
    ##---------------------------------------------------------------
    ##          Inputs for Splitting the Data on Infinity           -
    ##---------------------------------------------------------------



    def split_data_inputs(self):
        '''
        This function takes in inputs for how you want the split up the data
        so it can run faster on infinity. Generally, there are 5,000 galaxies
        in the NH simulation, so it is suggested to split the data evenly
        acoss nodes.
        
        '''
        
        parser = self.parser
        
        
        split_data_flags = {
            
            # This identifier will be used when naming the resulting output 
            # files:
            "--chunk_id":   {"type": int, "default": 0},
            
            # The number of chunks you want to split across nodes:
            "--n_chunks":   {"type": int, "default": 10},
            
    
            "--n_workers":  {"type": int, "default": 1},
            
            # How you want to split your data. You can split by galaxy ID or
            # filament ID
            "--split_type": {"type": str, "default": "gal_id"},
            
            
            # The number of k-nearest filament points that will be used to 
            # calculate a galaxy's distance from a filament:
            "--k_nearest_neighbors": {"type": int, "default": 20},
            
            # The minimum distance between interpolation points. These are 
            # expressed in Mpcs. The default is going to be set to 10 kpcs:
            "--min_interpolation_distance": {"type": float, "default": 10e-3}
            }
        
        # Adding all other paramenters into the parser:
        for flag, kwargs in split_data_flags.items():
            parser.add_argument(flag, **kwargs)
    
    
        
    
    ##---------------------------------------------------------------
    ##                      Print User Inputs                       -
    ##---------------------------------------------------------------


    
    
    def print_user_outputs(self):
        '''
        This function serves only to print out the inputs that were given.
        '''
        
        args = self.args
        
        print(f"The Time Step selected:\n{args.time_step}")
        print("*" * 50)
        print(f"The Bin Width selected:\n{args.bin_width}")
        print("*" * 50)
        print(f"Chunk ID selected:\n{args.chunk_id}")
        print(f"Number of chunks:\n{args.n_chunks}")
        print(f"Number of workers:\n{args.n_workers}")
        print("*" * 50)
        print(f"Type of splitting on the data:\n{args.split_type}")
        print("*" * 50)
        print(f"K-nearest neighbors to calculate a galaxy's distance from a" \
              f" filament:\n{args.k_nearest_neighbors}")
        print(f"The minimum interpolation distance:\n{args.min_interpolation_distance} Mpcs")
        
     
        
     
        
     

###########################################################################
###########################################################################
###                                                                     ###
###                FINDING THE DATA FOR THE CALCULATIONS                ###
###                                                                     ###
###########################################################################
###########################################################################



class RetrieveData:
    
    def __init__(self, data_path, filament_file_path, args):
        '''
        This class is responsible for retrieving the location of the galaxy and 
        filament data. 
        
        So far this function was created around the data found in the following
        file paths:
            
        Galaxy data: in TREE_STARS-Adapta_HOP
        
        Filament data: Filaments extracted from the gas, NOT the galaxies
        
        ----------------------------------------------------------------------
        
        Inputs:
            
            data_path (str):
                Taken galaxy data path taken in from 
            
            self.data_path (str):
                The data path where your galaxy data reside. This function was
                designed to intake galaxy data located inside the 
                TREE_STARS_Adapta_HOP file on Infinity. If running on a local
                machine, you can also utilize your own galaxy data, however
                ensure that your galaxies were pulled from the 
                TREE_STARS_Adapta_HOP file on Infinity.
                
            self.filament_file (str):
                The path where your gas filament file is located. 
            
        
        '''
        
        # Data path to where galaxies are located on Infinity:
        self.data_path = data_path
        
        # Galaxy or gas filament file location:
        self.filament_file = filament_file_path
        
        # Utilizing the user inputs from the UserInputs class:
        self.args = args
        
        
        # Initializing variables below that will be filled later depending on
        # the type of analysis:
        self.galaxy_list = None
        self.chunked_galaxies = None
        self.filament_critical_point_list = None
        self.filaments = None
        
        self.glob_galaxy_data()
        self.chunk_galaxies()
        self.load_filament_file()
    
    
    

    ##----------------------------------------------------------------
    ##                    Globbing the Galaxy Data                   -
    ##----------------------------------------------------------------



    def glob_galaxy_data(self):
        '''
        This function will glob up the galaxy data in the specified data path.
        It will return a list of the all galaxies found in that data path.
        '''
        self.galaxy_list = glob.glob(self.data_path + str("/gal_stars*"))
        
        
        
    
    ##----------------------------------------------------------------
    ##                  Reading in the Galaxy Data                   -
    ##----------------------------------------------------------------



    def chunk_galaxies(self):
        '''
        This function will split the data into n chunks, where n is the
        n_chunks parameter given by the user.
        '''
        
        chunks = self.args.n_chunks
        chunk_id = self.args.chunk_id
        
        self.galaxy_list.sort()
        galaxy_list = self.galaxy_list
        
        self.chunked_galaxies = np.array_split(galaxy_list, chunks)[chunk_id]
    
    
    ##----------------------------------------------------------------
    ##                    Reading in Filament File                   -
    ##----------------------------------------------------------------


    
    def load_filament_file(self):
        '''
        This function will read in either gas or galaxy filaments. The file
        must be in a NKDSkl file format. It will return a list of all the 
        node critical points found in the file and use those later for 
        linear interpolation and finding the galaxies that are closests to 
        these filaments.
        '''
        
        fil_dict = ft.get_fil_dict(self.filament_file)
        
        # Pulling out both filament information and critical point information:
        
        self.filaments = fil_dict["filaments"]
        
        self.filament_critical_point_list = fil_dict["critical_points"]
        

            
            
           
            
           
            
###########################################################################
###########################################################################
###                                                                     ###
###                    READ GALAXY AND FILAMENT DATA                    ###
###                                                                     ###
###########################################################################
###########################################################################


           
class ReadData:

    def __init__(self, 
                 galaxy_chunk_data, 
                 filament_data, 
                 critical_point_data, 
                 args):
        '''
        This function reads in the galaxy and filament data. This class will
        eat up most of your run time, especially reading in the galaxy data.
        
        Inputs:
            galaxy_chunk_data (str):
                The globbed up list of chunked galaxies to read in.
                
            filament_data (dict):
                The filament data location.
                
            critical_points_data (dict):
                
        '''
        
        self.imported_galaxy_chunk_data = galaxy_chunk_data
        self.imported_filament_data = filament_data
        self.imported_critical_point_data = critical_point_data
        
        self.args = args
        
        
        self.gal_stars_data = None
        self.galaxy_header_data = None
        self.filament_data = None
        self.filament_nodes = None
        self.combined_star_gal_data = None
        self.galaxy_data = None
        
        self.read_galaxy_data()
        self.create_filament_data()
        
        
            
    ##---------------------------------------------------------------
    ##        Reading Treebrick Information of the Galaxies         -
    ##---------------------------------------------------------------


    
    def read_galaxy_data(self):
        
        galaxy_list = self.imported_galaxy_chunk_data
        
        args = self.args


        self.galaxy_header_data = []
        self.gal_stars_data = []
        
        
        # Reading in the galaxies one at a time:
        for i in galaxy_list:
        
        # Because the galaxy data has the same naming convention, this is going
        # to find the number of the galaxy:
            match = re.search(r'_(\d+)$', i)
            number = int(match.group(1))
        
               # Actually read the galaxy data:
            obj = ReadGalData(gal_data = i,
                                alpha = args.alpha,
                                omega_m = args.omega_m,
                                omega_l = args.omega_l,
                                omega_k = args.omega_k,
                                aexp = args.aexp,
                                ntable = args.ntable,
                                h0 = args.H0)
        
            # Pull out the star level  information:
            star_info_df = obj.star_info
            
            star_info_df["Galaxy Number"] = number
            self.gal_stars_data.append(star_info_df)
            
        
            # Pull out the galaxy level information:
            galaxy_info_df = obj.galaxy_header_info
        
            # Append the galaxy information to gal_info
            self.galaxy_header_data.append(galaxy_info_df)
            
        
        self.combined_star_gal_data = pd.concat(self.gal_stars_data, ignore_index = True)

        self.combined_star_gal_data["Star_Mass (Solar Mass)"] = self.combined_star_gal_data["Star_Mass (Standard Mass Unit)"] * 1e11

        self.galaxy_data = pd.concat(self.galaxy_header_data, ignore_index = True)
        
        # Ensuring that galaxy number is an int64 instead of an object:
        self.galaxy_data = self.galaxy_data.astype({"Galaxy_Number": "int64"})



    
    ##---------------------------------------------------------------
    ##                Creating a Filament Data Frames               -
    ##---------------------------------------------------------------
    

    def create_filament_data(self):
        
        all_filaments = self.imported_filament_data
        critical_points = self.imported_critical_point_data
        
        filament_crt_points = []

        for i, fil in enumerate(all_filaments):
            
            cp1_idx = int(fil["cp1_idx"])
            cp2_idx = int(fil["cp2_idx"])
            
            ctr_points = pd.DataFrame({
                "Filament Index": [i],
                "cp1_idx": [cp1_idx],
                "cp2_idx": [cp2_idx],
                "cp1_type": [critical_points[cp1_idx]["cp_idx"]],
                "cp2_type": [critical_points[cp2_idx]["cp_idx"]], 
                "px1": [critical_points[cp1_idx]["px"]],
                "py1": [critical_points[cp1_idx]["py"]],
                "pz1": [critical_points[cp1_idx]["pz"]],
                "px2": [critical_points[cp2_idx]["px"]],
                "py2": [critical_points[cp2_idx]["py"]],
                "pz2": [critical_points[cp2_idx]["pz"]]
                
                })
            filament_crt_points.append(ctr_points)    
            
        self.filament_data = pd.concat(filament_crt_points, ignore_index = True)
        #self.filament_data = self.filament_data.query("(cp1_type == 3 or cp2_type == 3) and (cp1_type != 4 and cp2_type != 4)")
        
        
        
        #### EDIT: NOT SURE IF THIS CODE BELOW IS NEEDED, BUT I'M GOING TO
        #### KEEP IT HERE JUST IN CASE:
            
        # Next, combining together all of the node critical points:
        
        dd = defaultdict(list)

        for d in critical_points:
            for key, value in d.items():
                dd[key].append(value)

        # Now I'm going to transform dd into a Pandas dataframe:

        self.filament_nodes = pd.DataFrame.from_dict(dd).query("cp_idx == 3")
        
        

    
    
      
###########################################################################
###########################################################################
###                                                                     ###
###                 PERFORMING CALCULATIONS ON THE DATA                 ###
###                                                                     ###
###########################################################################
###########################################################################

  
class CalculateSFR:
    
    def __init__(self,
                 gal_stars_data,
                 galaxy_header_data,
                 args):
        
        self.star_data = gal_stars_data
        self.galaxy_header_data = galaxy_header_data      
        
        self.args = args
        
        
        self.age_of_universe = None
        self.combined_data = None
        
        
        # TESTING VARIABLES:
        self.testing_star_birth = None
        self.testing_adjusted_z = None
        self.testing_sfr_data = None
        self.testing_adjusted_sfr_data = None
        
        self.convert_time()
        self.combine_data_frames()
        
        
        
    ##---------------------------------------------------------------
    ##    Create a Function for Finding the Age of the Universe     -
    ##---------------------------------------------------------------


        
    def convert_time(self):
        '''
        Find the age of the Universe at the user-inputted time step
        '''
        
        args = self.args
        
        self.age_of_universe = CalculateMetrics.scale_factor_to_time(aexp = args.aexp,
                                                                     H0 = args.H0,
                                                                     omega_m = args.omega_m,
                                                                     omega_l = args.omega_l,
                                                                     omega_r = args.omega_r) # omega_r is negligable for this time step
        
        
            
    ##----------------------------------------------------------------
    ##                        Calculate SFHs                         -
    ##----------------------------------------------------------------

        
        
    def calculate_sfh(self):
        '''
        Calculate the star formation histories of the stars given the bins
        from the user input. This will calculate using both the star's age and
        the redshift adjusted star's age.
        '''
        
        print(f"AGE OF UNIVERSE:\n{self.age_of_universe}")
        args = self.args
        
        self.sfh_data = CalculateMetrics(data_frame = self.star_data,
                                         age_of_univ = self.age_of_universe).data_frame
        
        
        
    
        ##***************************************************************
        ##              Starting by Using the Star's Ages               *
        ##***************************************************************
    
    
        bin_width = args.bin_width
    
        rows = []
    
        for gal, filtered_galaxy_info in self.sfh_data.groupby("Galaxy Number"):
    
            bins = np.arange(0, filtered_galaxy_info["Star_Age (Years)"].max(), bin_width)
    
            # Count the number of stars in each bin:
            count_of_stars, _ = np.histogram(filtered_galaxy_info["Star_Age (Years)"], bins = bins)
    
            for edge, count in zip(bins[:-1], count_of_stars):
                rows.append({
                "Bins (Years)": edge,
                "Number of Stars in Bins": count,
                "Galaxy Number": gal
                })
    
    
        # Creating the data frame to use in plotting:
        star_birth = pd.DataFrame(rows)
        
        ########
        # TESTING CODE
        #######
        
        self.testing_star_birth = star_birth




        ##***************************************************************
        ##              Now Using the Adjusted Star's Age               *
        ##***************************************************************
    
    
        rows = []
    
        for gal, filtered_galaxy_info in self.sfh_data.groupby("Galaxy Number"): 
        
            # Star Age here is time of birth:
            bins = np.arange(0, filtered_galaxy_info["Lookback_Time_Adjusted_Z (Years)"].max(), bin_width)
        
            # Count the number of stars in each bin:
            count_of_stars, _ = np.histogram(filtered_galaxy_info["Lookback_Time_Adjusted_Z (Years)"], bins = bins)
        
            for edge, count in zip(bins[:-1], count_of_stars):
                rows.append({
                "Bins (Years)": edge,
                "Number of Stars in Bins": count,
                "Galaxy Number": gal
                })
            
        
        # Creating the data frame to use in plotting:
        star_birth_z_adjust = pd.DataFrame(rows)
        
        self.testing_adjusted_z = star_birth_z_adjust
        
        
        return star_birth, star_birth_z_adjust
    
    
    
    ##----------------------------------------------------------------
    ##                        Calculating SFRs                       -
    ##----------------------------------------------------------------


    def calculate_sfr(self):
        
        args = self.args
        
        sfh_data = self.sfh_data
        
        

        rows = []
    
        for gal, filtered_galaxy_info in sfh_data.groupby("Galaxy Number"):
    
            # Star Age here is time of birth:
            bins = np.arange(0, filtered_galaxy_info["Star_Age (Years)"].max(), args.bin_width)
    
            # Count the number of stars in each bin:
            count_of_stars, _ = np.histogram(filtered_galaxy_info["Star_Age (Years)"], bins = bins)
    
            # finding the mass per bin:
            mass_in_bin, _ = np.histogram(filtered_galaxy_info["Star_Age (Years)"],
                                      bins=bins,
                                      weights=filtered_galaxy_info["Star_Mass (Solar Mass)"])
    
    
            for edge, count, mass in zip(bins[:-1], count_of_stars, mass_in_bin):
                rows.append({
                "Bins (Years)": edge,
                "Number of Stars in Bins": count,
                "SFR (Solar Masses / Bin Width)": mass / args.bin_width,
                "Galaxy Number": gal
                })
    
    
        # Creating the data frame to use in plotting:
        sfr_data = pd.DataFrame(rows)
        
        self.testing_sfr_data = sfr_data
        
        



    
        ##***************************************************************
        ##              Now Using the Adjusted Star's Age               *
        ##***************************************************************
    
    
    
        rows = []
    
        for gal, filtered_galaxy_info in sfh_data.groupby("Galaxy Number"):
    
            # Creating bins of 50 Myrs:
            # Star Age here is time of birth:
            bins = np.arange(0, filtered_galaxy_info["Lookback_Time_Adjusted_Z (Years)"].max(), args.bin_width)
    
            # Count the number of stars in each bin:
            count_of_stars, _ = np.histogram(filtered_galaxy_info["Lookback_Time_Adjusted_Z (Years)"], bins = bins)
    
            # finding the mass per bin:
            mass_in_bin, _ = np.histogram(filtered_galaxy_info["Lookback_Time_Adjusted_Z (Years)"],
                                      bins=bins,
                                      weights=filtered_galaxy_info["Star_Mass (Solar Mass)"])
    
    
            for edge, count, mass in zip(bins[:-1], count_of_stars, mass_in_bin):
                rows.append({
                "Bins (Years)": edge,
                "Number of Stars in Bins": count,
                "SFR (Solar Masses / Bin Width)": mass / args.bin_width,
                "Galaxy Number": gal
                })
    
    
        # Creating the data frame to use in plotting:
        sfr_adjusted_data = pd.DataFrame(rows)
        
        self.testing_sfr_adjusted_data = sfr_adjusted_data
        
        return sfr_data, sfr_adjusted_data
    
    
    
        
    ##----------------------------------------------------------------
    ##                  Combining the SFR/SFH Data                   -
    ##----------------------------------------------------------------

    
    
    def combine_data_frames(self):
        
        star_birth, star_birth_z_adjusted = self.calculate_sfh()
        
        sfr_data, sfr_adjusted_data = self.calculate_sfr()
        
        
        # Renaming columns so when we join the tables together, there will be
        # no errors:
        sfh = star_birth[["Bins (Years)",
                          "Galaxy Number",
                          "Number of Stars in Bins"]].rename(
                              columns = {"Number of Stars in Bins": "SFH"}
                              )
        
        
        sfh_adjusted = star_birth_z_adjusted[["Bins (Years)",
                                              "Galaxy Number",
                                              "Number of Stars in Bins"]].rename(
                                                columns = {"Number of Stars in Bins": "SFH Z Adjusted"}
                                                )
        
        
        sfr = sfr_data[["Bins (Years)",
                        "Galaxy Number",
                        "SFR (Solar Masses / Bin Width)"]].rename(
                        columns = {"SFR (Solar Masses / Bin Width)": "SFR"}
                            )
        
        
        
        sfr_adjusted = sfr_adjusted_data[["Bins (Years)",
                                          "Galaxy Number",
                                          "SFR (Solar Masses / Bin Width)"]].rename(
                                              columns = {"SFR (Solar Masses / Bin Width)": "SFR Z Adjusted"}
                                              )
        
        
        combined_sfh = sfh.merge(sfh_adjusted,
                              on = ["Bins (Years)",
                                    "Galaxy Number"],
                              how = "outer")
        combined_sfr = combined_sfh.merge(sfr,
                                       on = ["Bins (Years)",
                                             "Galaxy Number"],
                                       how = "outer")
        
        self.combined_data = combined_sfr.merge(sfr_adjusted,
                                       on = ["Bins (Years)",
                                             "Galaxy Number"],
                                       how = "outer")
        
        self.combined_data["Bins (Gyrs)"] = self.combined_data["Bins (Years)"] / 1e9
        
        
        
        
        
        

   
############################################################################
############################################################################
###                                                                      ###
###                    CALCULATE DISTANCE TO FILAMENT                    ###
###                                                                      ###
############################################################################
############################################################################



class CalculateDistanceToFilament:
    
    def __init__(self,
                 filament_data,
                 critical_point_data,
                 galaxy_data,
                 sfr_data,
                 args):
        
        
        
        # Import data from previous classes:
        self.galaxy_data = galaxy_data
        self.sfr_data = sfr_data
        self.args = args
        self.filament_data = filament_data
        self.critical_point_data = critical_point_data
        
        
        
        # Set up data:
        self.filament_data_adjusted = None
        self.interpolated_points = None
        self.interpolated_filament_ids = None
        self.kdtree = None
        self.closest_distance_to_filament = None
        self.host_filament_id = None
        self.nearest_neighbors = None
        self.host_filament_critical_points = None
        
        
        
        # Final data after this class is ran:
        self.filament_galaxy_data = None
        
        
        
        # Run functions in the class:
        self.interpolate_between_critical_points()
        self.calculate_distance_to_filament()
        self.find_nearest_filament()
        self.combine_filament_and_gal_data()
        
        
       
    ##---------------------------------------------------------------
    ##                  Converting distances to Mpc                 -
    ##---------------------------------------------------------------

     
    def convert_distances(self, box_unit_distance):
        '''
        This function will convert the distance between points from standard
        box units to Mpcs. If you want to convert the entire box, you will
        have to add in a 0.5 offset!
        
        Check unit_l, its in centimeters. What you need to know is the size of
        your base simulation. 
        
        Look at the xmin and xmax for your gas cube! That will help you convert
        the distances.
        
        Try directly from the size of the box. What you currently have does not
        tell you the fraction of the gas cube you extracted. 
        
        The fraction is going to be useful later for analysis!
        '''
        
        unit_l = self.args.unit_l
        a_exp = self.args.aexp
        
        box_fraction = a_exp * self.args.boxlen * (self.args.xmax - self.args.xmin)
        
        # Convert the box units into Mpcs. This code was taken inside 
        # /data7b/NewHorizon/PP_NewH/create_listhalo.f90 while also finding 
        # the fraction size of my box compared to the NewHorizon box.
        
        return (box_unit_distance - 0.5) * box_fraction * unit_l / (1e6 * 3.08e18) 
        
        
        
        
    ##----------------------------------------------------------------
    ##                Interpolate Between the Points                 -
    ##----------------------------------------------------------------


    def interpolate_between_critical_points(self):  # rename to segments
        '''
        You have sampling points from disperse, use those! Maybe identify the 
        nearest 20 neighbors to then calculate the distanced a galaxy is. 
        You have to define the resolution. Is it worth
        '''
        
        filaments = self.filament_data
        args = self.args
        min_interpolation_distance = args.min_interpolation_distance
        
        
        filaments_units_adjusted = filaments.assign(px1 = self.convert_distances(filaments["px1"]),
                                                    py1 = self.convert_distances(filaments["py1"]),
                                                    pz1 = self.convert_distances(filaments["pz1"]),
                                                    px2 = self.convert_distances(filaments["px2"]),
                                                    py2 = self.convert_distances(filaments["py2"]),
                                                    pz2 = self.convert_distances(filaments["pz2"]))
        
        self.filament_data_adjusted = filaments_units_adjusted
            
        
        all_points = []
        all_filament_ids = []
        
        for _, rows in filaments_units_adjusted.iterrows():
            
            point_a = np.array([rows["px1"], rows["py1"], rows["pz1"]])
            point_b = np.array([rows["px2"], rows["py2"], rows["pz2"]])
            
            # First calculate the distance between sampling points in box units
            # then convert to Mpcs:
            
            distance_between_points= np.linalg.norm(point_b - point_a)
            
            
            #distance_between_points_Mpcs = self.convert_distances(distance_between_points_box_units)
            
            
            
            # Now check to see if the distance between the points is less than
            # 2 times the user-defined minimum. This is to ensure that the 
            # distance between each interpolated point is going to be at least
            # the length of the user inputted minimum:

            if distance_between_points < 2 * min_interpolation_distance:
                
                # If the distance between the sampling points are less than the
                # user defined minimum, don't interpolate. Just add the
                # sampling points to all_points and it's filament ID to
                # all_filament_ids:
                    
                all_points.extend([point_a, point_b])
                
                all_filament_ids.extend([rows["Filament Index"]] * 2)
            
            else:
                # If the distance between the sampling points is greater than
                # or equal to the user-inputted minimum distance, interpolate
                # points at the user defined minimum distance:
                
                number_of_interpolation_points = math.floor(distance_between_points / min_interpolation_distance)
                
                interpolate_between_crt_points = np.linspace(point_a, point_b, num = number_of_interpolation_points + 2) # Adding two here to account for the end points
                
                all_points.append(interpolate_between_crt_points)
                
                all_filament_ids.extend([rows["Filament Index"]] * len(interpolate_between_crt_points))
            
        self.interpolated_points = np.vstack(all_points)
        self.interpolated_filament_ids = np.array(all_filament_ids)



    
    ##----------------------------------------------------------------
    ##                Calculate Distance to Filament                 -
    ##----------------------------------------------------------------


    def calculate_distance_to_filament(self):
        '''
        This function will calculate the distance to a filament using KDTrees
        to speed up the runtime.
        '''
        
        
        self.kdtree = KDTree(self.interpolated_points)
        
        
    
    ##---------------------------------------------------------------
    ##                    Find closest filament                     -
    ##---------------------------------------------------------------



    def find_nearest_filament(self):
        '''
        Take KDTree results and find the closest filament based on the user-
        inputted k-nearest neighbors. This is to allow for a faster runtimes
        by calculating the distance a galaxy is away from the k-nearest
        neighbors instead of every single point in the simulation.
        '''
        
        galaxy_positions = self.galaxy_data[["Galaxy_Position_X (Mpc)", 
                                             "Galaxy_Position_Y (Mpc)", 
                                             "Galaxy_Position_Z (Mpc)"]].to_numpy()

        print("\nGalaxy position ranges:")
        print(self.galaxy_data[["Galaxy_Position_X (Mpc)","Galaxy_Position_Y (Mpc)","Galaxy_Position_Z (Mpc)"]].agg(["min","max"]))
        
                
        distances, nearest_index = self.kdtree.query(galaxy_positions, 
                                                     k = self.args.k_nearest_neighbors)
        
        # Now that we have the the nearest k-neighbors, we are now going to
        # find the the filament IDs and the actual closest distance a galaxy
        # is from a filament:
        
        potential_filaments = self.interpolated_filament_ids[nearest_index]
        
        filament_lookup_dictionary = {}
        
        for _, row in self.filament_data_adjusted.iterrows():
            point1 = np.array([row["px1"], row["py1"], row["pz1"]])
            point2 = np.array([row["px2"], row["py2"], row["pz2"]])
            
            filament_lookup_dictionary[row["Filament Index"]] = (point1, point2)
            
        
        # Creating empty arrays to fill in as we iterate through each galaxy
        # and calculate the distance to the filament, then finding the minimum
        # distance and assigning that filament as the galaxy's host:
            
        min_distance = np.empty(len(galaxy_positions))
        host_ids = np.empty(len(galaxy_positions), 
                                dtype = potential_filaments.dtype)
        
        for i, point in enumerate(galaxy_positions):
            
            # Some of the k-nearest neighbors can be apart of the same 
            # filament. This is to ensure that we are only looking at the 
            # unique filaments:
            
            unique_filaments = np.unique(potential_filaments[i])
            
            smallest_distance = np.inf
            best_filament = None
            
            for fil_id in unique_filaments:
                I, J = filament_lookup_dictionary[fil_id]
                
                real_distance = segment_distance(I, J, point)
                
                if real_distance < smallest_distance:
                    smallest_distance = real_distance
                    best_filament = fil_id
                
            
            min_distance[i] = smallest_distance
            host_ids[i] = best_filament
        
            
        # Convert to Mpcs for now... Not sure if I'll keep this or not.
        self.closest_distance_to_filament = min_distance
        
        self.host_filament_id = host_ids
        
        # Now I'm going to take the host_ids and find their node critical
        # points by merging the filament data onto the host filaments:
        
        host_critical_points = pd.DataFrame({
            "Filament Index": host_ids,
            })
        
        self.host_filament_critical_points = host_critical_points.merge(
            self.filament_data_adjusted,
            on = "Filament Index",
            how = "left"
            )
    
    
    
    ##---------------------------------------------------------------
    ##              Attach Host Filament to Galaxy Data             -
    ##---------------------------------------------------------------
    
    
    
    def combine_filament_and_gal_data(self):
        
        galaxy_number = self.galaxy_data["Galaxy_Number"].astype("int64")

        sfr_sfh_df = self.sfr_data.rename(columns = {"Galaxy Number": "Galaxy_Number"})                
        
        filament_results = pd.DataFrame({
            "Galaxy_Number": galaxy_number,
            "Distance To Filament (Mpcs)": self.closest_distance_to_filament,
            "Host Filament ID": self.host_filament_id,
            "cp1_type": self.host_filament_critical_points["cp1_type"],
            "cp1_idx": self.host_filament_critical_points["cp1_idx"],
            "cp2_type": self.host_filament_critical_points["cp2_type"],
            "cp2_idx": self.host_filament_critical_points["cp2_idx"],
            "cp1_px":  self.host_filament_critical_points["px1"],
            "cp1_py":  self.host_filament_critical_points["py1"],
            "cp1_pz":  self.host_filament_critical_points["pz1"],
            "cp2_px":  self.host_filament_critical_points["px2"],
            "cp2_py":  self.host_filament_critical_points["py2"],
            "cp2_pz":  self.host_filament_critical_points["pz2"]
            
            })
        
        
        galaxy_sfr_merge = self.galaxy_data.merge(sfr_sfh_df,
                                                  on = "Galaxy_Number",
                                                  how = "left")
        
        self.filament_galaxy_data = galaxy_sfr_merge.merge(filament_results,
                                                           on = "Galaxy_Number",
                                                           how = "left")
        
        


        
        
        
        
        
        
###########################################################################
###########################################################################
###                                                                     ###
###                      PRINT RESULTS TO INFINITY                      ###
###                                                                     ###
###########################################################################
###########################################################################


        
        
class PrintResults:
    
    def __init__(self,
                 args,
                 filament_galaxy_df):
        '''
        This function will print the data frame into a Parquet file. A 
        Parquet file will store the data in a smaller size, preventing an 
        overflow of storage when printing to Infinity (or your local machine).
        
        Inputs:
            self.args (parser):
                The user inputs given 
        '''
        
        self.args = args
        self.filament_galaxy_df = filament_galaxy_df
        
        self.convert_to_csv()
        
        
        
    def convert_to_csv(self):
        
        args = self.args
        
        self.filament_galaxy_df.to_parquet(
            f"{args.output_path}/sfh_sfr_ts{args.time_step}_bw{args.bin_width}_chunk{args.chunk_id}.zstd", 
            compression = "zstd")
        
        
        


############################################################################
############################################################################
###                                                                      ###
###                       CREATING A MAIN FUNCTION                       ###
###                                                                      ###
############################################################################
############################################################################


def main():
    
    # Start by getting the user inputs:
    user_inputs = UserInputs()
    args = user_inputs.args
    
    # Take those user inputs and retrieve the data locations:
    data_location = RetrieveData(data_path = args.data_path, 
                                 filament_file_path = args.filament_file_path, 
                                 args = args)
    
    # Now that we have the locations, read in the data into memory:
    read_data = ReadData(galaxy_chunk_data = data_location.chunked_galaxies, 
                         filament_data = data_location.filaments, 
                         critical_point_data = data_location.filament_critical_point_list, 
                         args = args)
    
    
    # Now that the data are in memory, calculate SFH and SFRs:
    sfr = CalculateSFR(gal_stars_data = read_data.combined_star_gal_data,
                       galaxy_header_data = read_data.galaxy_data,
                       args = args)
    
    # Calculate the distanced a galaxy is from the filament to determine its
    # host filametn:
    filament_distances = CalculateDistanceToFilament(filament_data = read_data.filament_data,
                                                     critical_point_data = data_location.filament_critical_point_list,
                                                     galaxy_data = read_data.galaxy_data,
                                                     sfr_data = sfr.combined_data, 
                                                     args = args)
    
    # Save the data to the disk:
    PrintResults(args = args,
                 filament_galaxy_df = filament_distances.filament_galaxy_data)

    
    

if __name__ == "__main__":
    main()
   
    
    
    
    
    


        
               
    
            
            
