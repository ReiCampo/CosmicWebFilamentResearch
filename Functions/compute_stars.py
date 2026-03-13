
############################################################################
############################################################################
###                                                                      ###
###                      READ GALAXY AND STARS DATA                      ###
###                                                                      ###
############################################################################
############################################################################

###  This function can read in the data found in                              
###  /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/GAL_XXXXX and     
###  separates the galaxy and star information in two different Pandas data   
###  frames.                                                                  

### If ever you need to check on the read_gal file that does unit conversion
### for you, you can look at read_gal_f90 in /data7b/NewHorizon/PP_NewH         


##----------------------------------------------------------------
##                  Importing Necessary Packages                 -
##----------------------------------------------------------------


import numpy as np
import pandas as pd
import time as t
from scipy.io import FortranFile


##----------------------------------------------------------------
##                        Begin Function                         -
##----------------------------------------------------------------


class ReadGalData:
    
    def __init__(self, gal_data, print_gal_variables = True):
        '''
        The initializing function for the class. This function will hold the
        header and galaxy/stars info in separate pandas data frames
        
        Inputs:
            self (self object):
                The object that will be called on within the ReadGalData class.
            
            gal_data: (New Horizons Data):
                The file path to your gal_stars_XXXXXX data. This function will
                take in the AdaptaHOP gal_stars_XXXXXX data found in the file
                path: 
                /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/GAL_XXXXX
                
            print_gal_varialbes (bool):
                The function will automatically be verbose and print out the 
                galaxy header information after the variables have been 
                declared. This may not be ideal if you are using this function
                to loop over multiple gal_stars_XXX files, so the function can 
                be silenced by setting this to False. 
                
        
        Outputs:
            self.header_info (Pandas data frame):
                The header information that contains the galaxy level info. See
                read_data function below to get more information on the columns
                contained in the data frame.
                
            self.star_info (Pandas data frame):
                The main data frame that contains the star information. For more
                details about the columns, see the read_data function below.
        '''
        
        self.gal_data = gal_data
        self.galaxy_header_info = pd.DataFrame()
        self.star_info = pd.DataFrame()
        self.verbose = print_gal_variables
        self.read_data()
        
    
    def read_data(self):
        '''
        This is the main function of ReadGalData. It will first read in the 
        Fortran binary file then initialize the galaxy variables in the
        appropriate order and save those to a pandas data frame. Next, it will
        read in the appropriate star information and then assign it to its own
        pandas data frame.
        
        The function is also verbose by printing out the galaxy header
        information.
        
        Inputs:
            self (self object):
                The object that will be called in the function to get the 
                self.gal_data, self.galaxy_header_info, and the self.star_info

        
        Outputs:
            self.galaxy_header_info (Pandas data frame):
                The galaxy information of your gal_stars_XXXXX data. The
                columns and their descriptions are as follows:
                
                Galaxy_Number (int):  The galaxy ID number. This ID can be
                used to "follow" the galaxy through time, however do not forget
                that mergers can happen so I'm fairly certain the ID number
                will change. (Will have to look into this later).
                
                Galaxy_Level (int):  The level the galaxy is. For example, it 
                can be classified as a halo, or sub halo, etc. As of writing
                this, I do not know what -1, 0, or 1 classify. My best guess
                is 1 is halo, 0 is subhalo and -1 I think is unidentified. (Will
                have to look into this later).
                
                Galaxy_Mass (Standard Unit):  The galaxy's mass in standard 
                units. You will have to multiply by 10^11 to get solar masses.
                
                Galaxy_Position_X (Mpc):  The X position of the galaxy.
                
                Galaxy_Position_Y (Mpc):  The Y position of the galaxy.
                
                Galaxy_Position_Z (Mpc):  The Z position of the galaxy.
                
                Galaxy_Velocity_X (km/s):  The velocity of the galaxy in the x
                direction.
                
                Galaxy_Velocity_Y (km/s):  The velocity of the galaxy in the y
                direction.
                
                Galaxy_Velocity_Z (km/s):  The velocity of the galaxy in the z
                direciton.
                
                Galaxy_Angular_Momentum_X (Standard Unit km s^-1):  The angular
                momentum of the galaxy in the x direction.
                
                Galaxy_Angular_Momentum_Y (Standard Unit km s^-1):  The angular
                momentum of the galaxy in the y direction.
                
                Galaxy_Angular_Momentum_Z (Standard Unit km s^-1):  The angular
                momentum of the galaxy in the Z direction.
                
                
            self.star_data (Pandas data frame):
                This data frame contains every star and their information in 
                the galaxy that you are looking at. The columns and their
                descriptions are as follows:
                
                Star_Position_X (Mpc):  The x position of the star. The units
                may be different here, will need to look into this later.
                
                Star_Position_Y (Mpc):  The y position of the star. The units
                may be different here, will need to look into this later.
                
                Star_Position_Z (Mpc):  The z position of the star. The units
                may be different here, will need to look into this later.
                
                Star_Velocity_X (km/s):  The velocity of the star in the x
                direction. Units may be different here, will have to check
                later.
                
                Star_Velocity_Y (km/s):  The velocity of the star in the y
                direction. Units may be different here, will have to check
                later.
                
                Star_Velocity_Z (km/s):  The velocity of the star in the z
                direction. Units may be different here, will have to check
                later.
                
                Star_Mass (Standard Mass Unit):  The mass of the star. I believe 
                that this still needs to be multiplied by *something* to get 
                solar masses, will need to look into this later.
                
                Star_Age (Standard Time Unit):  The age of the star. Will need 
                to look into how to change the units for this.
                
                Star_Metallicity (Dex):  The metallicity of the star. I am 
                assuming that the units are dex here, but will have to confirm
                this later.       

        '''
        
        # Setting up a timer to print later
        time1 = t.time()
        
        # Reading in the Fortran binary file:
        fortran_data = FortranFile(self.gal_data)
        
        
        # Reading in the variables in the same order found in the
        # compute_stars_kin9.f90:
        
        # These are values related to the galaxy. This will serve as my header
        # information. I am also squeezing the values since reading these in
        # will create a list, and I just want the values in those lists:
        my_number = np.squeeze(fortran_data.read_record("i"))
        level = np.squeeze(fortran_data.read_record("i"))
        m = np.squeeze(fortran_data.read_record("d"))
        px, py, pz = np.squeeze(fortran_data.read_record("d", "d", "d"))
        vx, vy, vz = np.squeeze(fortran_data.read_record("d", "d", "d"))
        Lxg, Lyg, Lzg = np.squeeze(fortran_data.read_record("d", "d", "d"))
        
        # Number of stars in the galaxy:
        npart = np.squeeze(fortran_data.read_record("i"))
        
        if self.verbose == True:
        
            print(f"Number of Stars in Galaxy {my_number}: {npart}")
            
            print("-" * 50)
   
            print(f"Galaxy Number:\n{my_number}\n" \
                f"Galaxy Level:\n{level}\n" \
                f"Mass of Galaxy (Standard Units): {m}\n" \
                f"Galaxy Position X (Mpc):\n{px}\n" \
                f"Galaxy Position Y (Mpc):\n{py}\n" \
                f"Galaxy Position Z (Mpc):\n{pz}\n" \
                f"Galaxy Velocity X (km/s):\n{vx}\n" \
                f"Galaxy Velocity Y (km/s):\n{vy}\n" \
                f"Galaxy Velocity Z (km/s):\n{vz}\n" \
                f"Galaxy Angular Momentum X (Standard Unit km s^-1):\n{Lxg}\n" \
                f"Galaxy Angular Momentum Y (Standard Unit km s^-1):\n{Lyg}\n" \
                f"Galaxy Angular Momentum Z (Standard Unit km s^-1):\n{Lzg}\n")
            
            print("*" * 50)
            
            
        self.galaxy_header_info = pd.DataFrame({
            
            "Galaxy_Number": [my_number],
            "Galaxy_Level": [level],
            "Galaxy_Mass (Standard Unit)": [m],
            "Galaxy_Position_X (Mpc)": [px],
            "Galaxy_Position_Y (Mpc)": [py],
            "Galaxy_Position_Z (Mpc)": [pz],
            "Galaxy_Velocity_X (km/s)": [vx],
            "Galaxy_Velocity_Y (km/s)": [vy],
            "Galaxy_Velocity_Z (km/s)": [vz],
            "Galaxy_Angular_Momentum_X (Standard Unit km s^-1)": [Lxg],
            "Galaxy_Angular_Momentum_Y (Standard Unit km s^-1)": [Lyg],
            "Galaxy_Angular_Momentum_Z (Standard Unit km s^-1)": [Lzg]   
            
        })
        
        # These are values pertaining to the stars.
        x_stars = fortran_data.read_record("f")
        y_stars = fortran_data.read_record("f")
        z_stars = fortran_data.read_record("f")
        
        vx_stars = fortran_data.read_record("f")
        vy_stars = fortran_data.read_record("f")
        vz_stars = fortran_data.read_record("f")
        
        mass_stars = fortran_data.read_record("f")
        ids = fortran_data.read_record("f") # The star ids are being thrown away in the original Fortran file, so I will assign for now but will not use it later
        age_stars = fortran_data.read_record("f")
        zz_stars = fortran_data.read_record("f")
        
        
        self.star_info = pd.DataFrame({
            
            # Star Information Within the Single Galaxy:
            "Stars_Position_X (Mpc)": x_stars,
            "Stars_Position_Y (Mpc)": y_stars,
            "Stars_Position_Z (Mpc)": z_stars,
            "Star_Velocity_X (km/s)": vx_stars,
            "Star_Velocity_Y (km/s)": vy_stars,
            "Star_Velocity_Z (km/s)": vz_stars,
            "Star_Mass (Standard Mass Unit)": mass_stars,
            "Star_Age (Standard Time Unit)": age_stars,
            "Star_Metallicity (Dex)": zz_stars
        })
        
        time2 = t.time()
        
        print("*" * 50)
        print(f"Data Read Time: {(time2 - time1)} seconds")
        print("*" * 50)
        
        return self.galaxy_header_info, self.star_info
    
    