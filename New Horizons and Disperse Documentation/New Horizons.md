### What is New Horizons?

New Horizons is a hydrodynamical cosmological simulation that can create high resolution sims without sacrificing statistical power. This is done by taking a large box, simulating it, then zooming in about 16 Mpc^3 and re-simulating that region. This allows for greater resolution and can capture the interstellar medium with some statistics.

New Horizons can model gas cooling, star formation, feedback from massive stars, and feedback from supermassive black holes.

### New Horizons Data for Cosmic Web Group


New Horizons data is stored on the Infinity supercomputer which the group uses to run simulations on. The files are in the format of binary files (I'm pretty sure this is right?), so Fortran and Python scripts have been written to read in the data. 

If you are creating a script to read in any New Horizons data, variable declaration matters! Let's say the file contains:

### New Horizons Meeting Notes:

Data sets:
New Horizons
Horizon AGN

New Horizon: Snapshot from redshit 50 to 0, well sampled in between, like at 6, 7, and 0. Not as many snapshots as timesteps, this is for memory purpose. You have to select which snapshots you want. Tehre are more than 1000 snapshots, most of them are below redshift 5. There is a good sample here but depends on the timestep, whats a lower redshift is better than at higher reshift. There are just as many snapshots at z < 1 as z > 1. Bias for lower redshift. The smallest redshift is z = 0.17. 

Most of our analysis is between redshift 0 and 2, so the sampling is good!

In one snapshot, you have RAMSES files.
1 snapshot is one directory, and contained are different files corresponding to different components. 
There is one file per processer to run the simulation. The sim was ran with more than 5000 cores. For each type of file, you have 5000 cores! Each processor spits out its own file.

### Raw files

So in one directory you have particle files, contains the particles which can be the dark matter, the stars, and the black holes. For these you ahve at least m, position, v, for stars you have date of birht, metalicity, etfc.  These are simple to use b/c they are just files

In NH, we have specific files for black holes (we may or may not use these)

-----------------------------------
Then you have the hydro files (the properties as the gas). This is not computed as particles, but as a grid in an adaptive fashion. Think back to the grid back in your earlier notes. You can view this in memory by using a tree. An *AMR* structure is a tree where that's the root node of the tree and splits into quadrants. You store all levels, and the values at any given level, they are the average of the values below (see picture). It's a octree, at any given level the pixels all have the same volume. This is unlike KDTrees. You want this in simulation b/c you want to be as able to simulate any part!

You have position, identity of parent node, and identity of child nodes. The top is called the root (again see picture). This is a complex structure to store in memory!! We have files just to re cord the structure of the nodes!! No hydrodynamics has been applied. AMR files just contain the structure of the tree.

The hydrodynamic files contain the structure of the variables, density, temperature, pressure, metallicity, and velocity of the gas in the cell. (for each node)

grav files - This stores the gravitational potential in each node.
sink files - stores black hole information. This is given by every processor

You will have to stitch all of this together if you want all the stars within the radius of a halo, so you'll have to stitch together all the nodes that are related to what you're lookign at.

Uses a hilbert curve to split the job evenly between the 5000 processors. At every timestep we try to even out hte jobs. The hilbert curve covers the plane, a 1D structure that covers a 2D area if you push them to infinity. It refines through a multi-scale process. The length of hte curve scales with the stuff you're looking at. The length will scale with the amount of work you have to do. You'll split it up into 5000 pieces and redo it every time step. This is a purely analytical. Each processor stores a hilbert key, the function is in a fortran routine to retreive the certain process within an area/volume. So it's nice to only have to look in one space instead of the ENTIRE file. 

In practice, we rarely use the raw files. But if we want to look at intra-halo light you have to go back to the raw files and do calculations. 

There is one single file called the info file, recalculates the base information. Cosmological parameters, scale factor, size of hte box, age of the universe, etc. The number attached to the info_XXX.txt file is trivial, has no meaning. You can use these file to check if thigns are in the redshift range you want. User readable! All other files mentioned above are fortran binary files.

particle files - you have routines that use particle and hydrodynamic files. 


### Postprocessed files

These have their own assumptions built into them. 

treebricks - Galaxy or halo catalogs. For halo you use dm particles, for galaxy you use star particles. They are all produced with a structure finder called HaloMaker. Each of these files have their own little tricks that we add as we go. However there are version that just focus on halos and ones that focus on galaxies. HaloMaker is a clustering algorithm that works on density. It tries to find particles that cluster together, it uses a "friend of friend " method to find this. (fof). The idea is you are going to link the particles whether they are close or not. You start with any particle and link it to particles taht are close enough by choosing a kernal that is fixed or determined by local properties of hte field. There are thresholds (so if there are 4 particles in a cluster, that gets rejected). The standard is rejecting any cluster < 50. THIS IS USER DEFINED. 20% of the linking length is a tuned parameter that is used in the kernal for fof algorithms. 

AdaptaHOP is a method to find substructures. So it can find sub halos within halos. WE used tweed et all 2009 to do this. It's reasonably efficient to find structures, but it is not good for galaxy mergers. So it's pretty bad at identifying the mass of each component in a merger. Charlotte developed a method called Velociraptor. If you want to identify structure, you wil lhave to use density and in phase space
So when you see 6DFOF, that means it's a friend of friend in pahse space. 

To link galaxies through time, you will have to link the treebricks. So you can use TreeMaker. whe n you track a halo back in time, you track it back to all its progenitors. You track back the particles. If two particles have been lost somewhere else, you are goign to trash those particles that left. Out of all the projenitors, we consider the most massive one is the main one, the previous halo. 

Tree bricks are available already. They have been heavily tuned and they are good. But the actual tree that philogony of haloes is not something that is readily available, you have to run treemaker. You give treemaker that gives you a snapshot of the link and it creates a geneology tree. The output structure is a tree, not easy to store in Fortran, but we have base routines to read them. 

So you ahve separate tools of your halos and galaxies and a tool to help you link the galaxies you want through time. 

Since treebricks is a fortran binary, it is a list (a catalog), so for each galaxy/halo you have the identity, level (structure, substructure, etc), host (like next sub halo, 0 or -1 if they don't make sense), mass, positions, velocities, spin, shape, size, etc... The more fancy metrics are not necessarily calculated the way you want them to be. But they are good for first checks. But generally you want to calculate your own spin, shape, size, etc... Treebricks are in mega parsecs, km/s, mass is in standard unit, you'll have to multiply by 10^11 for solar masses. 

If you want to extrac t filaments from the galaxies, all you need are the treebricks. 

Hoewver we are missing the list and properties of particles! We need these particles for things like star formation rates!

The new version of HaloMaker, it spits out the particles for hte list of galaxies outomatically. This is not in the treebrick files!!

### GAL_XXX Directories


Gal is a directory of files, one directory per snapshot. 

gal_XXX_IDGAL.dat, and you have one file per galaxy. For each galaxy in this file, you have the list of the stars. It is also a fortran binary catalog where you have galaxy IDs, mass, but then you have a list of stars, positions, velocities, age, metallickty, etc... 

So you can start with a treebrick to find the identities of galaxies with a certain mass, then you can read in gal_XXX data and the stars in those galaxies. So you can use this information to do ground truth analysis, or a mock studies. If you wanted to start with mock studies, you will take the list of stars, or dm particles if you're just looking at the halo, but dm particles aremnt the most informative. You can refine stallar properties based on these files. Sometimes the stars captured in tehse files are actually background stars and aren't actually in teh galaxy! 

tbirth is stored in weird code units. 





To compute sfh's, we can do it in the "ground truth" way. sort stars by age and reconstruct the sfr from this as a function of time. You can use mocks by passing it through a filter, so each star particle is assigned a spectrum and passed through a dust screen then you can use a filter (like LSST or HST) and it can give you u-band luminosity. Dust is a big degeneracy factor, so you may want to do something more proper. You should use SKiRT, this takes into account scattering in all directions. You can cast millions of rays through your gas, they scatter with certain probabilities, etc... You need external code when you have ray scatter, absorption is never really a problem. You can have resonance lines that add a bit more energy. SKiRT only uses dust, not resonance. NH is may not be resolved enough for RASCAS yet. 

We don't follow dust in cosmological simulations. WE first have to understand the ground truth, get the filaments from the galaxies. Think about how you're going to build your routine to build the sfh. 

For RAMSES, it is a mesh code, SPH codes fair super bad with shocks, which isn't good for galaxy evolution! 

output_dir is where the raw data files are stored

### Info_XXX.txt

levelmin- lowest refinement, 
levelmax - highest refinement
nstep_corase - corase time step

boxln - in code units, default 1
time - in weird code units
aexp - expansion factor
h0 - scale factor
unit l - length
unit d - density
unit t - time, can use this stuff to conversions

order type - hilbert curve used to  decompose work
underneath this, thesea re the hilbert keys. You need the function to read this 

tree_dm - halo treebricks
tree_stars_adapta - This is the reference one we use. 
tree_stars - Lower refinement

For each tree brick, you have a corresponding gal directory

gal_stars_IDofGALAXIES, 5000 or so galaxie are in NH. 

# Today:

We focused on treebricks and gal_stars_XXX file. 

Next: running python scripts on HPC queue
Using treemaker 
Reconstructing filaments with disperse
Extracting cubes of gas


3/19/2026 Research Meeting Notes
There are about 5000 galaxies in New Horizons

To calculate the distance from a filament:
Remember Disperse gives you samples along the saddle point. in between the points you have segments, which have a center (this is just by taking the average between the extremeties). You can do distance to center of segment, that's reasonable enough to find the distance from the filament. Sometimes, though, you can get galaxies that are a part of one filament, but when you calculate teh distance from the center of the segment, sometimes you may get that the galaxy is closer to a segment in another filament than the one that the galaxy actually is in. You can take the vector between the two points of the segment and you can determine if the position of hte galaxy is negative. If it is, the closest point is A (see image). What we really need is the normal vector. Once you have the normal vector, you can get the ex, ey, and ez of the vector. To get the orthogonal distance, you take the dot product between two vectors and see if it is equal to zero. Zero means orthogonality. You will have to select which vector is equal to 1 and equal to 0. Then you have to solve for the third component of the vector. Once you have those two vectors, to get the third vector, you take the cross product of those. Then you project AG onto the e2 e3 plane and calculate the length of the resultant vector. 

KDTree is super fast at calculating the distance from points. It avoid computing all of hte point. It first segements the the cell wehre you are and searches the cell that you are in (and maybe the cells near it). This only works when you have point to point. 