import os
from os.path import join
import re
import numpy as np
import time
import struct
import matplotlib.pyplot as plt
import pandas as pd
from scipy.io import FortranFile
import ctypes as c

class ReadHaloShapes:
    def __init__(self,file_path=None):
        
        """
        Make a treebricks dictionary out of file
        """
        self.file_path = file_path
        self.HaloShapes_dict = None
        self.grid = None
        self.read_data()

    def read_data(self):
        t0 = time.time()
        f = FortranFile(self.file_path, 'r')

        nbhaloes,nn = f.read_record('i')
        print('nhaloes:',nbhaloes)
        nshells = int((nn-3)/16)
        print('nshells:',nshells)

        self.HaloShape_dict = {}
        self.HaloShape_dict['nbhaloes'] = nbhaloes
        self.HaloShape_dict['nshells'] = nshells
        self.HaloShape_dict['haloes'] = []

        grid = f.read_record('f4').reshape(nbhaloes,nn, order='F')

          
        for i in range(nbhaloes):

            halo = grid[i]

            #each halo gets its own dictionary
            halo_dict = {}

            id = int(halo[0])
            mtot = halo[1]
            rvir = halo[2]

            halo_dict['ID'] = id
            halo_dict['mtot'] = mtot
            halo_dict['rvir'] = rvir

            for j in range(nshells):
                #for each shell, make a dictionary for shell vals
                halo_dict[f'shell {j+1}'] = []

            for j in range(nshells):


                shell_dict = {}

                a = halo[3+j*16+0]
                b = halo[3+j*16+1]
                c = halo[3+j*16+2]

                Lx = halo[3+j*16+3]
                Ly = halo[3+j*16+4]
                Lz = halo[3+j*16+5]

                aa_x = halo[3+j*16+6]
                aa_y = halo[3+j*16+7]
                aa_z = halo[3+j*16+8]

                bb_x = halo[3+j*16+9]
                bb_y = halo[3+j*16+10]
                bb_z = halo[3+j*16+11]

                cc_x = halo[3+j*16+12]
                cc_y = halo[3+j*16+13]
                cc_z = halo[3+j*16+14]

                m = halo[3+j*16+15]


                shell_dict['a'] = a
                shell_dict['b'] = b
                shell_dict['c'] = c
                shell_dict['Lx'] = Lx
                shell_dict['Ly'] = Ly
                shell_dict['Lz'] = Lz
                shell_dict['aa_x'] = aa_x
                shell_dict['aa_y'] = aa_y
                shell_dict['aa_z'] = aa_z
                shell_dict['bb_x'] = bb_x
                shell_dict['bb_y'] = bb_y
                shell_dict['bb_z'] = bb_z
                shell_dict['cc_x'] = cc_x
                shell_dict['cc_y'] = cc_y
                shell_dict['cc_z'] = cc_z
                shell_dict['m'] = m

                halo_dict[f'shell {j+1}'].append(shell_dict)

            self.HaloShape_dict['haloes'].append(halo_dict)




        #HaloShape_dict --> halo_dict --> shell dict

        return self.HaloShape_dict

class ReadTreebrick_highp:
    #low precision Treebricks (HAGN)
    def __init__(self,file_path=None):
        
        """
        Make a treebricks dictionary out of file
        """
        self.file_path = file_path
        self.treebricks_dict = None
        
        #***# RLC Edit:
        # Adding in a separate data frame because I feel that the 
        # treebricks_dict is not actually that useful. I will keep 
        # treebricks_dict stored in its own variable, but halo_data will 
        # contain the actual values needed to do analysis and create plots:
        self.halo_data = pd.DataFrame()
        self.read_data()
   
        
    def read_data(self):


        t0 = time.time()
        f = FortranFile(self.file_path, 'r')
        
        n_bodies = f.read_record('i')
        m_part = f.read_record('f')
        a_exp = f.read_record('f')
        om_t = f.read_record('f')
        age_universe = f.read_record('f')
        n_sub = f.read_record('i')
        n_haloes = n_sub[0] + n_sub[1]
        
        
        #***# RLC Edit: Starting off by creating a pandas dataframe instead of a 
        #***# dictionary:
        
        self.treebricks_dict = pd.DataFrame()
        
        # Mutating in the columns:
        self.treebricks_dict = self.treebricks_dict.assign(
            nbodies = n_bodies,
            mpart = m_part,
            aexp = a_exp,
            omega_t = om_t,
            age_univ = age_universe,
            nhost = n_sub[0],
            nsub = n_sub[1],
            nhaloes = n_haloes)
        
        
        # RLC Edit: Commenting this out for now before I get a working version
        # of this class going:
        
        # self.treebricks_dict = {}
        # self.treebricks_dict['nbodies'] = nbodies
        # self.treebricks_dict['mpart'] = mpart
        # self.treebricks_dict['aexp'] = aexp
        # self.treebricks_dict['omega_t'] = omega_t
        # self.treebricks_dict['age_univ'] = age_univ
        # self.treebricks_dict['nhost'] = nsub[0]
        # self.treebricks_dict['nsub'] = nsub[1]
        # nhaloes = nsub[0]+nsub[1]
        # self.treebricks_dict['nhaloes'] = nhaloes
        # #initialize empty list to hold all halos
        # self.treebricks_dict['haloes'] = []
        
        print('nbodies:', n_bodies, 
              '\nmpart:', m_part,
              '\naexp:', a_exp,
              '\nomega_t:', om_t,
              '\nage:', age_universe,
              '\nnsub:', n_sub,
              '\nnhaloes:', n_haloes)
        
        
        #***# RLC Edit: This variable is not used in the rest of the class so I
        # am going to comment it out for now:
        # nb_of_haloes = self.treebricks_dict['nhost'] + self.treebricks_dict['nsub']
        
        #define the length of the box in Mpc for New Horizon
        lbox_nh = 2.0 * 0.07 * (142.8571428) * a_exp 
        lbox_hagn = (142.8571428) * a_exp
        halfbox = lbox_nh / 2.0
        
        # Adding the box length in to the data frame:
        self.treebricks_dict = self.treebricks_dict.assign(
            lbox_NH = lbox_nh,
            lbox_HAGN = lbox_hagn
        )
        
        # Commenting this out for now in case I need it later:
        # self.treebricks_dict['lbox_NH'] = lbox_NH
        # self.treebricks_dict['lbox_HAGN'] = lbox_HAGN
        
        
        #def read_halo():
        t1 = time.time()
        for i in range(n_haloes):
            halo_dict = {}

            npart = f.read_record('i')
            halo_dict['npart'] = npart.tolist()

            #members = f.read_record('i').reshape(npart, order='F')
            members = f.read_record('i')
            halo_dict['members'] = members.tolist()

            my_number = f.read_record('i')
            halo_dict['my_number'] = my_number.tolist()

            my_timestep = f.read_record('i')
            halo_dict['my_timestep'] = my_timestep.tolist()


            level_ids = f.read_record('i')
            level_ids = level_ids.tolist()

            level, host_halo, host_sub, nchild, nextsub = level_ids[0], level_ids[1], level_ids[2], level_ids[3], level_ids[4]

            halo_dict['level'] = level
            halo_dict['host_halo']  = host_halo
            halo_dict['host_sub'] = host_sub
            halo_dict['nchild'] = nchild
            halo_dict['nextsub'] = nextsub

            mass = f.read_record('f') #d for NH
            # RLC Edit: Going to get rid of mass.tolist() just to see what happens
            halo_dict["mass"] = mass
            #halo_dict['mass'] = mass.tolist()

            p = f.read_record('f')
            p = p.tolist()
            
            # RLC Edit: I think this line below is an error? It starts with 
            # py = p[0], but based on the other variables like vx, vx = v[0] not
            # vy = v[0]...  For now, I am going to say px = p[0] then ask about 
            # this later.]
            # px and py are swapped in some catlogs. Revert back to the old way
            # Create a check to see if there's a match between 
            # py,px,pz = p[0],p[1],p[2]
            py, px, pz = p[0], p[1], p[2]
            halo_dict['px'] = px
            halo_dict['py'] = py
            halo_dict['pz'] = pz

            v = f.read_record('f')
            v = v.tolist()
            vx,vy,vz = v[0],v[1],v[2]  
            halo_dict['vx'] = vx
            halo_dict['vy'] = vy
            halo_dict['vz'] = vz

            L = f.read_record('f')
            L = L.tolist()
            Lx,Ly,Lz = L[0],L[1],L[2]
            halo_dict['Lx'] = Lx
            halo_dict['Ly'] = Ly
            halo_dict['Lz'] = Lz

            shape = f.read_record('f')
            shape = shape.tolist()
            rmax,a,b,c = shape[0],shape[1],shape[2],shape[3]
            halo_dict['rmax'] = rmax
            #write old shapes as well
            halo_dict['a'] = a
            halo_dict['b'] = b
            halo_dict['c'] = c
            
            energy = f.read_record('f')
            energy = energy.tolist()
            ek,ep,et = energy[0],energy[1],energy[2]
            halo_dict['ek'] = ek
            halo_dict['ep'] = ep
            halo_dict['et'] = et

            spin = f.read_record('f')
            halo_dict['spin'] = spin.tolist()

            virial = f.read_record('f')
            virial = virial.tolist()
            rvir,mvir,tvir,cvel = virial[0],virial[1],virial[2],virial[3]
            halo_dict['rvir'] = rvir
            halo_dict['mvir'] = mvir
            halo_dict['tvir'] = tvir
            halo_dict['cvel'] = cvel 

            halo_profile = f.read_record('f')
            halo_profile = halo_profile.tolist()
            rho_0, r_c = halo_profile[0],halo_profile[1]
            halo_dict['rho_0'] = rho_0
            halo_dict['r_c'] = r_c
            
            #Positions are in Mpc, we now put them back in "code units", that is assuming the length of the simulation is 1.
            #for New Horizon scale
            """
            
            px_NH = (halo_dict['px']/self.treebricks_dict['lbox_NH']) + 0.5
            halo_dict['px_NH'] = px_NH
            
            py_NH = (halo_dict['py']/self.treebricks_dict['lbox_NH']) + 0.5
            halo_dict['py_NH'] = py_NH
           
            pz_NH = (halo_dict['pz']/self.treebricks_dict['lbox_NH']) + 0.5
            halo_dict['pz_NH'] = pz_NH
            
            rvir_NH = (halo_dict['rvir']/self.treebricks_dict['lbox_NH'])
            halo_dict['rvir_NH'] = rvir_NH

            
            #for Horizon AGN scale
            px_HAGN = (halo_dict['px']/self.treebricks_dict['lbox_HAGN']) + 0.5
            halo_dict['px_HAGN'] = px_HAGN
            
            py_HAGN = (halo_dict['py']/self.treebricks_dict['lbox_HAGN']) + 0.5
            halo_dict['py_HAGN'] = py_HAGN
            
            pz_HAGN = (halo_dict['pz']/self.treebricks_dict['lbox_HAGN']) + 0.5
            halo_dict['pz_HAGN'] = pz_HAGN
            
            rvir_HAGN = (halo_dict['rvir']/self.treebricks_dict['lbox_HAGN'])
            halo_dict['rvir_HAGN'] = rvir_HAGN
            """

            
            #return halo_dict
            self.halo_data = self.halo_data._append(halo_dict, ignore_index = True)
            
            
        #t1 = time.time()
        #for i in range(nhaloes):
            #halo_dict = read_halo()
            #self.treebricks_dict['haloes'].append(halo_dict)
            
        t2 = time.time()
        print('Reading haloes took {:0.2f} secs.'.format(t2-t1))
        print('Total time was {:0.2f} secs.'.format(t2-t0))
            
        return self.halo_data, self.treebricks_dict
    
class ReadTreebrick_lowp:
    #low precision Treebricks (to read halos)
    def __init__(self, file_path = None, haloshape_dict = None):
        
        """
        Make a treebricks dictionary out of file
        """
        self.file_path = file_path
        self.haloshape_dict = haloshape_dict
        self.treebricks_dict = None
        
        #***# RLC Edit: Just like in ReadTreebrick_highp, I will create a
        # separate halo dictionary. I will retain the treebrick_dict, but that
        # will not be the primary data frame I will use to analyze and plot the
        # haloes.
        self.halo_data = pd.DataFrame()
        self.read_data()
   
        
    def read_data(self):


        def get_halo_shape(halo_dict,ID):
            nbhaloes = halo_dict['nbhaloes']
            ids = [halo_dict['haloes'][i]['ID'] for i in range(nbhaloes)]
            ids = np.asarray(ids)
    
            halo_idx = np.where(ids == ID)
            idx = int(halo_idx[0][0])

    
            a = halo_dict['haloes'][idx]['shell 3'][0]['a']
            b = halo_dict['haloes'][idx]['shell 3'][0]['b']
            c = halo_dict['haloes'][idx]['shell 3'][0]['c']
    
    
            return a,b,c


        t0 = time.time()
        f = FortranFile(self.file_path, 'r')
        
        n_bodies = f.read_record('i')
        m_part = f.read_record('f')
        a_exp = f.read_record('f')
        om_t = f.read_record('f')
        age_universe = f.read_record('f')
        n_sub = f.read_record('i')
        n_haloes = n_sub[0] + n_sub[1]
        
    
        self.treebricks_dict = pd.DataFrame()
        
        self.treebricks_dict = self.treebricks_dict.assign(
            nbodies = n_bodies,
            mpart = m_part[0],
            aexp = a_exp[0],
            omega_t = om_t[0],
            age_univ = age_universe[0],
            nh_old = n_sub[0],
            nsub_old = n_sub[1],
            nhaloes = n_haloes
        )
    
        
        print('nbodies:', n_bodies,
              '\nmpart:', m_part,
              '\naexp:', a_exp,
              '\nomega_t:', om_t,
              '\nage:', age_universe,
              '\nnsub:',n_sub,
              '\nnhaloes:', n_haloes)
        
        print("=" * 50)
        print("Passed Creating a Dictionary")
        print("=" * 50)
        
        nb_of_haloes = self.treebricks_dict['nh_old'] + self.treebricks_dict['nsub_old']
        
        #define the length of the box in Mpc for New Horizon
        lbox_nh = 2.0 * 0.07 * (142.8571428) * a_exp[0] 
        lbox_hagn = (142.8571428) * a_exp[0]
        halfbox = lbox_nh / 2.0
        
        # Assigning the new columns to the treebricks_dict dataframe:
        self.treebricks_dict = self.treebricks_dict.assign(
            lbox_NH = lbox_nh,
            lbox_HAGN = lbox_hagn
        )
        
        #**# RLC Edit: Commenting this block out for now in case I need it 
        # later:
        # self.treebricks_dict['lbox_NH'] = lbox_NH
        # self.treebricks_dict['lbox_HAGN'] = lbox_HAGN
        
        
        #def read_halo():
        t1 = time.time()
        for i in range(n_haloes):
            halo_dict = {}

            npart = f.read_record('i')
            halo_dict['npart'] = npart.tolist()

            #members = f.read_record('i').reshape(npart, order='F')
            members = f.read_record('i')
            halo_dict['members'] = members.tolist()

            my_number = f.read_record('i')
            halo_dict['my_number'] = my_number.tolist()

            my_timestep = f.read_record('i')
            halo_dict['my_timestep'] = my_timestep.tolist()


            level_ids = f.read_record('i')
            level_ids = level_ids.tolist()
            print(f"Level IDS: {level_ids}")

            level,host_halo,host_sub,nchild,nextsub = level_ids[0],level_ids[1],level_ids[2],level_ids[3],level_ids[4]
            halo_dict['level'] = level
            halo_dict['host_halo']  = host_halo
            halo_dict['host_sub'] = host_sub
            halo_dict['nchild'] = nchild
            halo_dict['nextsub'] = nextsub

            mass = f.read_record('f') #d for NH
            halo_dict['mass'] = mass.tolist()

            p = f.read_record('f')
            p = p.tolist()
            py,px,pz = p[0],p[1],p[2]
            halo_dict['px'] = px
            halo_dict['py'] = py
            halo_dict['pz'] = pz

            v = f.read_record('f')
            v = v.tolist()
            vx,vy,vz = v[0],v[1],v[2]  
            halo_dict['vx'] = vx
            halo_dict['vy'] = vy
            halo_dict['vz'] = vz

            L = f.read_record('f')
            L = L.tolist()
            Lx,Ly,Lz = L[0],L[1],L[2]
            halo_dict['Lx'] = Lx
            halo_dict['Ly'] = Ly
            halo_dict['Lz'] = Lz

            shape = f.read_record('f')
            shape = shape.tolist()
            rmax,a_,b_,c_ = shape[0],shape[1],shape[2],shape[3]
            halo_dict['rmax'] = rmax
            #write old shapes as well
            halo_dict['old_a'] = a_
            halo_dict['old_b'] = b_
            halo_dict['old_c'] = c_

            ### OVERWRITE HAGN haloshapes with the shapes calculated by Charlotte and replace them 
            ### -------------------------------------------------------------
            ### RLC Edit: Make sure to come back to this to find the correct
            ### haloshape_dict in order for this section to work!
            # a,b,c = get_halo_shape(self.haloshape_dict,ID = my_number.tolist())
            # halo_dict['a'] = a 
            # halo_dict['b'] = b
            # halo_dict['c'] = c

            energy = f.read_record('f')
            energy = energy.tolist()
            ek,ep,et = energy[0],energy[1],energy[2]
            halo_dict['ek'] = ek
            halo_dict['ep'] = ep
            halo_dict['et'] = et

            spin = f.read_record('f')
            halo_dict['spin'] = spin.tolist()

            virial = f.read_record('f')
            virial = virial.tolist()
            rvir,mvir,tvir,cvel = virial[0],virial[1],virial[2],virial[3]
            halo_dict['rvir'] = rvir
            halo_dict['mvir'] = mvir
            halo_dict['tvir'] = tvir
            halo_dict['cvel'] = cvel 

            halo_profile = f.read_record('f')
            halo_profile = halo_profile.tolist()
            rho_0, r_c = halo_profile[0],halo_profile[1]
            halo_dict['rho_0'] = rho_0
            halo_dict['r_c'] = r_c
            

            #add halo dict to main dict under key of my_number
            
            #return halo_dict
            self.halo_data = self.halo_data._append(halo_dict, ignore_index = True)
            
            
        t2 = time.time()
        print('Reading haloes took {:0.2f} secs.'.format(t2-t1))
        print('Total time was {:0.2f} secs.'.format(t2-t0))
            
        return self.halo_data, self.treebricks_dict
    
class GalaxyCatalog:
    #High precision Treebricks to read galaxies 
    def __init__(self,file_path=None):
        
        """
        Make a treebricks dictionary out of file
        """
        self.file_path = file_path
        self.treebricks_dict = pd.DataFrame()
        self.gal_data = pd.DataFrame()
        self.read_data()
   
        
    def read_data(self):
        t0 = time.time()
        f = FortranFile(self.file_path, 'r')
        
        n_bodies = f.read_record('i')
        m_part = f.read_record('d')
        a_exp = f.read_record('d')
        om_t = f.read_record('d')
        age_universe = f.read_record('d')
        n_sub = f.read_record('i')
        n_max = n_sub[0] + n_sub[1]
        
        self.treebricks_dict = self.treebricks_dict.assign(
            nbodies = n_bodies,
            mpart = m_part,
            aexp = a_exp,
            omega_t = om_t,
            age_univ = age_universe,
            nhost = n_sub[0],
            nsub = n_sub[1],
            nmax = n_max)
        
        # self.treebricks_dict = {}
        # self.treebricks_dict['nbodies'] = nbodies
        # self.treebricks_dict['mpart'] = mpart
        # self.treebricks_dict['aexp'] = aexp
        # self.treebricks_dict['omega_t'] = omega_t
        # self.treebricks_dict['age_univ'] = age_univ
        # self.treebricks_dict['nb_of_galaxies'] = nsub[0]
        # self.treebricks_dict['nb_of_subgals'] = nsub[1]
        # nmax = nsub[0]+nsub[1]
        # self.treebricks_dict['nmax'] = nmax
        # #initialize empty list to hold all halos
        # self.treebricks_dict['galaxies'] = []
        
        print('nbodies:',n_bodies,
              'mpart:', m_part,
              'aexp:', a_exp,
              'omega_t:', om_t,
              'age:', age_universe,
              'nsub:',n_sub,
              'nmax:', n_max)
        
        
        
                
        #define the length of the box in Mpc for New Horizon
        lbox_nh = 2.0 * 0.07 * (142.8571428) * a_exp 
        lbox_hagn = (142.8571428) * a_exp
        halfbox = lbox_nh / 2.0
        
        self.treebricks_dict = self.treebricks_dict.assign(
          lbox_NH = lbox_nh,
          lbox_HAGN = lbox_hagn  
        )
        # self.treebricks_dict['lbox_NH'] = lbox_NH
        # self.treebricks_dict['lbox_HAGN'] = lbox_HAGN
        
        
        #def read_halo():
        t1 = time.time()
        for i in range(n_max):
            gal_dict = {}

            npart = f.read_record('i')
            gal_dict['npart'] = npart.tolist()

            #members = f.read_record('i').reshape(npart, order='F')
            members = f.read_record('i')
            gal_dict['members'] = members.tolist()

            my_number = f.read_record('i')
            gal_dict['my_number'] = my_number.tolist()

            my_timestep = f.read_record('i')
            gal_dict['my_timestep'] = my_timestep.tolist()


            level_ids = f.read_record('i')
            level_ids = level_ids.tolist()
            #print('level ids',level_ids)
            level,host_gal,host_subgal,nchild,nextsub = level_ids[0],level_ids[1],level_ids[2],level_ids[3],level_ids[4]
            gal_dict['level'] = level
            gal_dict['host_gal']  = host_gal
            gal_dict['host_subgal'] = host_subgal
            gal_dict['nchild'] = nchild
            gal_dict['nextsub'] = nextsub

            mass = f.read_record('d') #d for NH
            #print('mass',mass)
            gal_dict['mass'] = mass.tolist()

            p = f.read_record('d')
            #print('p',p)
            p = p.tolist()
            py,px,pz = p[0],p[1],p[2]
            gal_dict['px'] = px
            gal_dict['py'] = py
            gal_dict['pz'] = pz


            v = f.read_record('d')
            #print('v',v)
            v = v.tolist()
            vx,vy,vz = v[0],v[1],v[2]  
            gal_dict['vx'] = vx
            gal_dict['vy'] = vy
            gal_dict['vz'] = vz

            L = f.read_record('d')
            #print('L',L)
            L = L.tolist()
            Lx,Ly,Lz = L[0],L[1],L[2]
            gal_dict['Lx'] = Lx
            gal_dict['Ly'] = Ly
            gal_dict['Lz'] = Lz

            shape = f.read_record('d')
            #print('shape',shape)
            shape = shape.tolist()
            rmax,a,b,c = shape[0],shape[1],shape[2],shape[3]
            gal_dict['rmax'] = rmax
            gal_dict['a'] = a 
            gal_dict['b'] = b
            gal_dict['c'] = c

            energy = f.read_record('d')
            #print('energy',energy)
            energy = energy.tolist()
            ek,ep,et = energy[0],energy[1],energy[2]
            gal_dict['ek'] = ek
            gal_dict['ep'] = ep
            gal_dict['et'] = et

            spin = f.read_record('d')
            #print('spin',spin)
            gal_dict['spin'] = spin.tolist()

            sigma = f.read_record('d')
            #print('sigma',sigma)
            sigma = sigma.tolist()
            sig,sigma_bulge,m_bulge = sigma[0],sigma[1],sigma[2]
            gal_dict['sigma'] = sig
            gal_dict['sigma_bulge'] = sigma_bulge
            gal_dict['m_bulge'] = m_bulge
            
            virial = f.read_record('d')
            virial = virial.tolist()
            rvir,mvir,tvir,cvel = virial[0],virial[1],virial[2],virial[3]
            gal_dict['rvir'] = rvir
            gal_dict['mvir'] = mvir
            gal_dict['tvir'] = tvir
            gal_dict['cvel'] = cvel 

            
            halo_profile = f.read_record('d')
            #print('halo profile',halo_profile)
            halo_profile = halo_profile.tolist()
            rho_0, r_c = halo_profile[0],halo_profile[1]
            gal_dict['rho_0'] = rho_0
            gal_dict['r_c'] = r_c
            
            #stellar density profiles
            
            nbins = f.read_record('i')
            #print('nbins:',nbins)
            
            rr = f.read_record('d') #array (nbins)
            gal_dict['rr'] = rr
            rho = f.read_record('d')
            gal_dict['rr'] = rho
            
            #print('shape of density profile',np.shape(rho))
            


            #add gal dict to main dict 
            
            #return halo_dict
            self.gal_data = self.gal_data._append(gal_dict, ignore_index = True)
            
            
        t2 = time.time()
        print('Reading galaxies took {:0.2f} secs.'.format(t2-t1))
        print('Total time was {:0.2f} secs.'.format(t2-t0))
            
        return self.gal_data, self.treebricks_dict