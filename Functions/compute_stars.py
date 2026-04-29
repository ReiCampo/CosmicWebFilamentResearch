
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

   
###  The catalogs to view all the galaxy information is in                        
###  /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/tree_bricksXXX.       
###  Make sure to use the same time step as the Gal file! So for example if you   
###  pulled a galaxy from GAL_00970, make sure your tree bricks file says         
###  tree_bricks970!                                                              
                                                 

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
    
    def __init__(self, 
                 gal_data, 
                 alpha,
                 omega_m,
                 omega_l,
                 omega_k,
                 aexp,
                 ntable,
                 h0,
                 print_gal_variables = True):
        '''
        The initializing function for the class. This function will hold the
        header and galaxy/stars info in separate pandas data frames
        
        Inputs:
            self (self object):
                The object that will be called on within the ReadGalData class.
            
            gal_data (New Horizons Data):
                The file path to your gal_stars_XXXXXX data. This function will
                take in the AdaptaHOP gal_stars_XXXXXX data found in the file
                path: 
                /data7b/NewHorizon/TREE_STARS_AdaptaHOP_dp_SCnew_gross/GAL_XXXXX
                
            alpha (int):
                Alpha will be used when calling convert_time() later in the 
                code. This is the integration accuracy (or integration step
                size) when calculating the conformal time derivative of the
                scale factor and the cosmic time derivative of the scale 
                factor. In the original Fortran code, it was hard coded to be
                1e-6.
                
            omega_m (double):
                The Omega_matter cosmological paramter found in the
                info_XXXXX.txt files (in the file path /data7b/NewHorizon/INFO 
                on infinity). You should be using the omega_m found in your 
                specific New Horizons snapshot info_XXXXX.txt file!
                
            omega_l (double):
                The Omega_lambda cosmological paramter found in the
                info_XXXXX.txt files (in the file path /data7b/NewHorizon/INFO 
                on infinity). You should be using the omega_l found in your 
                specific New Horizons snapshot info_XXXXX.txt file!
                
            omega_k (double):
                The Omega_k (curvature) cosmological paramter found in the
                info_XXXXX.txt files (in the file path /data7b/NewHorizon/INFO 
                on infinity). You should be using the omega_k found in your 
                specific New Horizons snapshot info_XXXXX.txt file!
                
            aexp (double):
                The cosmological scale factor found in the info_XXXXX.txt files 
                (in the file path /data7b/NewHorizon/INFO on infinity). You 
                should be using the aexp found in your specific New Horizons 
                snapshot info_XXXXX.txt file!
                
            ntable (int):
                The resolution of the lookup table. This determines how many
                subsamples of the Friedmann equations to take, then create a 
                fixed-size table with evenly spaced entries. This was hard-
                coded in the Fortran file to be 1000.
                
            h0 (double):
                The Hubble constant corresponding to your timestep. Make sure
                to use the Hubble constant found in your info_XXXXX.txt file!
                
                
            print_gal_varialbes (bool):
                The function will automatically be verbose and print out the 
                galaxy header information after the variables have been 
                declared. This may not be ideal if you are using this function
                to loop over multiple gal_stars_XXX files, so the function can 
                be silenced by setting this to False. 
                
        
        Outputs:
            self.header_info (Pandas data frame):
                The header information that contains the galaxy level info. See
                read_data() function below to get more information on the columns
                contained in this data frame.
                
            self.star_info (Pandas data frame):
                The main data frame that contains the star information. For more
                details about the columns, see the read_data() function below.
        '''
        
        self.gal_data = gal_data
        self.galaxy_header_info = pd.DataFrame()
        self.star_info = pd.DataFrame()
        self.verbose = print_gal_variables
        self.alpha = alpha
        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_k = omega_k
        self.aexp = aexp
        self.ntable = ntable
        self.h0 = h0
        
        self.read_data()
        self.convert_time()
        
    
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
                may be different here, will need to look into this later. These
                values have not been recentered onto the parent galaxy.
                
                Star_Position_Y (Mpc):  The y position of the star. The units
                may be different here, will need to look into this later. These
                values have not been recentered onto the parent galaxy.
                
                Star_Position_Z (Mpc):  The z position of the star. The units
                may be different here, will need to look into this later. These
                values have not been recentered onto the parent galaxy.
                
                Star_Velocity_X (km/s):  The velocity of the star in the x
                direction. Units may be different here, will have to check
                later. These values have not been recentered onto the parent
                galaxy.
                
                Star_Velocity_Y (km/s):  The velocity of the star in the y
                direction. Units may be different here, will have to check
                later. These values have not been recentered onto the parent
                galaxy.
                
                Star_Velocity_Z (km/s):  The velocity of the star in the z
                direction. Units may be different here, will have to check
                later. These values have not been recentered onto the parent
                galaxy.
                
                Star_Mass (Standard Mass Unit):  The mass of the star. I believe 
                that this still needs to be multiplied by *something* to get 
                solar masses, will need to look into this later.
                
                Star_Age (Standard Time Unit):  The age of the star in code
                units. Can be converted to years with the function 
                convert_time().
                
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
        x_stars = fortran_data.read_record("d")
        y_stars = fortran_data.read_record("d")
        z_stars = fortran_data.read_record("d")
        
        vx_stars = fortran_data.read_record("d")
        vy_stars = fortran_data.read_record("d")
        vz_stars = fortran_data.read_record("d")
        
        mass_stars = fortran_data.read_record("d")
        ids = fortran_data.read_record("i") # The star ids are being thrown away in the original Fortran file, so I will assign for now but will not use it later
        age_stars = fortran_data.read_record("d")
        zz_stars = fortran_data.read_record("d")
        
        
        self.star_info = pd.DataFrame({
            
            # Star Information Within the Single Galaxy:
            "Star_Position_X (Mpc)": x_stars,
            "Star_Position_Y (Mpc)": y_stars,
            "Star_Position_Z (Mpc)": z_stars,
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
    
    
    def convert_time(self):
        '''
        This function will convert the time unit of the star's age into years.
        This function was originally written in Fortran and was translated
        using Claude.ai
        
        Inputs:
            self (self object):
                The self object that will be called throughout the function.
                
        Outputs:
            self.star_info (Pandas data frame):
                This will return a Pandas data frame with all of the above 
                columns found in the read_data() function, but with an added
                column:
                
                Star_Age (Years): The age of the stars in years. This was 
                calculated using the conversion functions found in the 
                Fortran code. 
        
        '''
        
        
        def dadtau(axp, O_mat_0, O_vac_0, O_k_0):
            """da/dtau — conformal time derivative of scale factor"""
            return axp**2 * np.sqrt(O_mat_0/axp**3 + O_vac_0 + O_k_0/axp**2)

        def dadt(axp, O_mat_0, O_vac_0, O_k_0):
            """da/dt — cosmic time derivative of scale factor"""
            return axp * np.sqrt(O_mat_0/axp**3 + O_vac_0 + O_k_0/axp**2)

        def friedman(O_mat_0, O_vac_0, O_k_0, alpha, axp_min, ntable):
            """
            Reproduce the Fortran friedman subroutine.
            
            Parameters:
                O_mat_0 : float - matter density parameter (Omega_m)
                O_vac_0 : float - vacuum/dark energy density (Omega_lambda)
                O_k_0   : float - curvature density (Omega_k), often 0
                alpha   : float - step size accuracy (e.g. 1e-4)
                axp_min : float - minimum scale factor to integrate to (e.g. 1e-4)
                ntable  : int   - number of entries in output lookup table
            
            Returns:
                axp_out, hexp_out, tau_out, t_out — each of length ntable+1
            """

            # --- Pass 1: count total steps to determine nskip ---
            axp_tau = 1.0
            axp_t   = 1.0
            tau = 0.0
            t   = 0.0
            nstep = 0

            while axp_tau >= axp_min or axp_t >= axp_min:
                nstep += 1

                dtau = alpha * axp_tau / dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0)
                axp_tau_pre = axp_tau - dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0) * dtau / 2.0
                axp_tau = axp_tau - dadtau(axp_tau_pre, O_mat_0, O_vac_0, O_k_0) * dtau
                tau -= dtau

                dt = alpha * axp_t / dadt(axp_t, O_mat_0, O_vac_0, O_k_0)
                axp_t_pre = axp_t - dadt(axp_t, O_mat_0, O_vac_0, O_k_0) * dt / 2.0
                axp_t = axp_t - dadt(axp_t_pre, O_mat_0, O_vac_0, O_k_0) * dt
                t -= dt

            age_tot = -t
            nskip = nstep // ntable
            print(f"Age of the Universe (in unit of 1/H0) = {age_tot:.3e}")

            # --- Pass 2: build the lookup table ---
            axp_out  = np.zeros(ntable + 1)
            hexp_out = np.zeros(ntable + 1)
            tau_out  = np.zeros(ntable + 1)
            t_out    = np.zeros(ntable + 1)

            axp_tau = 1.0;  tau = 0.0
            axp_t   = 1.0;  t   = 0.0
            nstep = 0;      nout = 0

            t_out[0]    = t
            tau_out[0]  = tau
            axp_out[0]  = axp_tau
            hexp_out[0] = dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0) / axp_tau

            while axp_tau >= axp_min or axp_t >= axp_min:
                nstep += 1

                dtau = alpha * axp_tau / dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0)
                axp_tau_pre = axp_tau - dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0) * dtau / 2.0
                axp_tau = axp_tau - dadtau(axp_tau_pre, O_mat_0, O_vac_0, O_k_0) * dtau
                tau -= dtau

                dt = alpha * axp_t / dadt(axp_t, O_mat_0, O_vac_0, O_k_0)
                axp_t_pre = axp_t - dadt(axp_t, O_mat_0, O_vac_0, O_k_0) * dt / 2.0
                axp_t = axp_t - dadt(axp_t_pre, O_mat_0, O_vac_0, O_k_0) * dt
                t -= dt

                if nstep % nskip == 0:
                    nout += 1
                    if nout <= ntable:
                        t_out[nout]    = t
                        tau_out[nout]  = tau
                        axp_out[nout]  = axp_tau
                        hexp_out[nout] = dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0) / axp_tau

            # Force last entry
            t_out[ntable]    = t
            tau_out[ntable]  = tau
            axp_out[ntable]  = axp_tau
            hexp_out[ntable] = dadtau(axp_tau, O_mat_0, O_vac_0, O_k_0) / axp_tau

            return axp_out, hexp_out, tau_out, t_out
        
        # --- Build the lookup tables ---
        aexp_out, _, tau_frw, t_frw = friedman(O_mat_0 = self.omega_m, 
                                        O_vac_0 = self.omega_l, 
                                        O_k_0 = self.omega_k, 
                                        ntable = self.ntable,
                                        alpha = 1e-6, 
                                        axp_min = 1e-3)
        
        time_simu = np.interp(self.aexp, aexp_out[::-1], tau_frw[::-1])
        
        # --- Interpolate all star ages at once ---
        age_conformal = self.star_info["Star_Age (Standard Time Unit)"].values
        time_cosmic = np.interp(age_conformal, tau_frw[::-1], t_frw[::-1])
        
        # --- Convert to physical age in years and store ---
        self.star_info["Star_Age (Years)"] = (
            (time_simu - time_cosmic) / (self.h0 * 1e5 / 3.08e24) / (365.0 * 24.0 * 3600.0)
        )
        
        return self.star_info

    
class RecenterStars:
    

    def __init__(self, gal_info, star_info):
        '''
        This class recenters the stars to the host galaxy's center of mass.
        This is because when you first use ReadGalData, the stars' positions
        and velocities are not centered on the host galaxy's center of mass, so
        this class should be used immediately after ReadGalData.
        
        Inputs:
            self (self object): 
                The object that will be called within the RecenterStars class.
                
            gal_info (Pandas Data Frame):
                The galaxy header information of the host galaxy. This should
                be the same galaxy header information returned from 
                ReadGalData.
            
            star_info (Pandas Data Frame):
                The star information of the host galaxy. This should be the 
                same star info returned from ReadGalData.

        '''
        
        self.gal_info = gal_info
        self.star_info = star_info
        self.adjusted_star_info = star_info
        self.recenter_stars()
        
    
    def recenter_stars(self):
        '''
        This is the main function of RecenterStars. This will take in the
        star's positions and velocities and adjust them based on the host
        galaxy's positions and velocities.
        
        Inputs:
            self (self object):
                The object that will be called inside the recenter_stars 
                function
                
        Outputs:
            adjusted_star_info (Pandas Data Frame):
                The star info data frame that contains the recentered star
                positions and velocities. The following adjusted columns are
                as follows:
                    
                    Star_Position_X (Mpc):  Subtracts the host galaxy's x
                    position from the star's x position.
                    
                    Star_Position_Y (Mpc):  Subtracts the host galaxy's y
                    position from the star's y position.
                    
                    Star_Position_Z (Mpc):  Subtracts the host galaxy's z
                    position from the star's z position.
                    
                    Star_Velocity_X (km/s):  Subtracts the host galaxy's x
                    velocity from the star's x velocity.
                    
                    Star_Velocity_Y (km/s):  Subtracts the host galaxy's y
                    velocity from the star's y velocity.
                    
                    Star_Velocity_Z (km/s):  Subtracts the host galaxy's z
                    velocity from the star's z velocity

        '''
        
        # Only finding the galaxy_id_number to use later in the print 
        # statement:
        galaxy_id_number = self.gal_info["Galaxy_Number"].values[0]
        
        # Now finding the galaxy's com positions and velocities:
        galaxy_com_positions = [self.gal_info["Galaxy_Position_X (Mpc)"].values[0],
                                 self.gal_info["Galaxy_Position_Y (Mpc)"].values[0],
                                 self.gal_info["Galaxy_Position_Z (Mpc)"].values[0]]
        
        galaxy_com_velocities = [self.gal_info["Galaxy_Velocity_X (km/s)"].values[0],
                                 self.gal_info["Galaxy_Velocity_Y (km/s)"].values[0],
                                 self.gal_info["Galaxy_Velocity_Z (km/s)"].values[0]]
        
        
        
        print(f"Galaxy {galaxy_id_number} center of mass coordinates in Mpc:\n" \
              f"{galaxy_com_positions}")
        print("*" * 50)
        print(f"Galaxy {galaxy_id_number} center of mass velocities in km/s:\n" \
              f"{galaxy_com_velocities}")
        print("*" * 50)
  
    
        # This is the main section of the function. This recenters the star's 
        # positions and velocities:
        self.adjusted_star_info = self.adjusted_star_info.assign(**{
                                                                    "Star_Position_X (Mpc)": self.star_info["Star_Position_X (Mpc)"] - galaxy_com_positions[0],
                                                                    "Star_Position_Y (Mpc)": self.star_info["Star_Position_Y (Mpc)"] - galaxy_com_positions[1],
                                                                    "Star_Position_Z (Mpc)": self.star_info["Star_Position_Z (Mpc)"] - galaxy_com_positions[2],
                                                                    "Star_Velocity_X (km/s)": self.star_info["Star_Velocity_X (km/s)"] - galaxy_com_velocities[0],
                                                                    "Star_Velocity_Y (km/s)": self.star_info["Star_Velocity_Y (km/s)"] - galaxy_com_velocities[1],
                                                                    "Star_Velocity_Z (km/s)": self.star_info["Star_Velocity_Z (km/s)"] - galaxy_com_velocities[2]})
        
        return self.adjusted_star_info
    
    

class CalculateLineOfSight:
    
    def __init__(self, gal_info, adjusted_star_info):
        
        self.gal_info = gal_info
        self.star_info = adjusted_star_info
        self.rotated_positions = None
        self.rotated_velocities = None
        self.los_calculation()
        
    def los_calculation(self):
        
        star_positions = np.column_stack([self.star_info["Star_Position_X (Mpc)"].values,
                          self.star_info["Star_Position_Y (Mpc)"].values,
                          self.star_info["Star_Position_Z (Mpc)"].values])
        
        star_velocities = np.column_stack([self.star_info["Star_Velocity_X (km/s)"].values,
                           self.star_info["Star_Velocity_Y (km/s)"].values,
                           self.star_info["Star_Velocity_Z (km/s)"].values])
        
        star_mass = self.star_info["Star_Mass (Standard Mass Unit)"].values
        
        
        # Finding the distance each star is from the galaxy's center:
        
        galaxy_com_positions = [self.gal_info["Galaxy_Position_X (Mpc)"].values[0],
                                 self.gal_info["Galaxy_Position_Y (Mpc)"].values[0],
                                 self.gal_info["Galaxy_Position_Z (Mpc)"].values[0]]
        
        distance_to_center = np.column_stack([star_positions[:, 0] - galaxy_com_positions[0],
                                           star_positions[:, 1] - galaxy_com_positions[1],
                                           star_positions[:, 2] - galaxy_com_positions[2]])
        
        
        # Calculating the angular momentum of each star from the center of the
        # galaxy:
            
        L_total = np.sum(star_mass[:, None] * np.cross(star_positions, star_velocities), axis = 0)
        
        # Normalizing so I can get the direction of the vector normal to the 
        # angular momentum:
            
        L_normal = L_total / np.linalg.norm(L_total)
        
        
        # Now I need to create the new x and y axis, so I will do that by using
        # an arbitrary vector like we discussed in research meeting:
            
        helper_vector = [1, 0, 0]
        
        # Creating the x axis by crossing the helper vector and the L_normal
        x_cross = np.cross(helper_vector, L_normal)
        x_axis = x_cross / np.linalg.norm(x_cross)
        
        # Creating the y axis by crossing the x_axis with the L_normal:
        y_axis = np.cross(L_normal, x_axis)
        
        # Now that I have my new x, y, and z axis, I will now create the
        # rotation matrix.
        
        rotation_matrix = np.array([x_axis, y_axis, L_normal])
        
        # Now let's actually rotate the position of the stars:
        
        self.rotated_positions = pd.DataFrame({
                    "Rotated_Position_X (Mpc)": (rotation_matrix @ star_positions.T).T[:, 0],
                    "Rotated_Position_Y (Mpc)": (rotation_matrix @ star_positions.T).T[:, 1],
                    "Rotated_Position_Z (Mpc)": (rotation_matrix @ star_positions.T).T[:, 2]
                    })

        self.rotated_velocities = pd.DataFrame({
                    "Rotated_Velocity_X (km/s)": (rotation_matrix @ star_velocities.T).T[:, 0],
                    "Rotated_Velocity_Y (km/s)": (rotation_matrix @ star_velocities.T).T[:, 1],
                    "Rotated_Velocity_Z (km/s)": (rotation_matrix @ star_velocities.T).T[:, 2]
                    })
        
        return self.rotated_positions, self.rotated_velocities
    

        
        
        
        