# Multi-objective evolutionary algorithm for optimising the geometry of a top-hinged wave energy converter
These are the codes that were used to run simulations for Edwards et al, 'The effect of device geometry on the performance of a wave energy converter'

MOEA_control.py
-	This file is how the multi-objective optimisations are run
-	In this file, you specify the possible values of the geometric parameters you want to consider
-	It calls the ‘MOEA.py’ function to find the population and Pareto Front (which it plots at the end of the code)
-	It also calls ‘plot_2D.py’ to plot the side-view 2D plots of the shapes on the Pareto Front
  
MOEA.py
-	This is the multi-objective evolutionary algorithm. See MOEA_control for the input variables. The outputs are final population and Pareto Front (vectors, F5 and kW for each)
-	It defines an initial population (randomly, within the limits of each parameter). For each shape in the initial population, the shape is tested to make sure the curves do not cross over (which would be of course infeasible) using ‘test_shapes’
-	To find the objective functions (F5 and kW), ‘given_geom_vector_find_kW_F_and_I55’ is called. If I55 is negative, the vector is thrown away
-	It finds the Pareto Front by extracting nondominated individuals from the initial population. To do so, it calls ‘extract_nondom_vectors’
-	Then, the algorithm loop is started
o	One parent is chosen from current population, using ‘select_parent_1’
o	The other parent is chosen from the current PF randomly
o	The child is created using the parents, using ‘create_offspring’
o	Child is added to the population and, if it is nondominated, it is added to the PF
o	This is repeated for the specified number of generations

plot_2D.py
-	This file turns the geometry vector into a plot of the 2D shape (x-z)
test_shapes.py
-	This looks at a theoretical shape based on the vector to see if it is feasible (if it crosses, it is not feasible)

given_geom_vector_find_kW_F_and_I55.py
-	First, the center of buoyancy of the shape is found using ‘calculate_xCB_and_zCB,’ the surface area is found using ‘calculate_surface_area’ and volume is calculated using ‘calculate_volume’
-	Then, the hydrodynamic coefficients are found by running WAMIT using ‘WAMIT_wraper_func’
-	Then, the relevant parameters are calculated (F5, kW, I55)

extract_nondom_vectors.py
-	This code forms the PF from the population. The inputs are vectors, parameter 1 values for the population and parameter 2 values for the population
-	It keeps track of dominated solutions (to then output nondominated solutions)
-	To determine if a solution is dominated, it uses ‘dominance_btwn_two_individuals’

select_parent_1.py
-	This follows the typical way of selecting a parent for a MOEA (running a ‘tournament’)

create_offspring.py
-	This follows the typical way of creating an offspring for a MOEA, using crossover and mutation
-	Once a vector is defined, the shape is tested (using ‘test_shapes’) and, if it’s okay, the relevant parameters are found (using ‘given_geom_vector_find_kW_F_and_I55’)

calculate_xCB_and_zCB.py
-	This find the x and z coordinates of the center of the shape (which will be the center of buoyancy)

calculate_surface_area.py
-	This finds surface area of the WEC

calculate_volume.py
-	This find the volume of the WEC

WAMIT_wrapper_func.py
-	This is the main file that controls running WAMIT.
-	There is an option to run a lookup_table, so that you can use what you’ve already run
-	Assuming you need to run WAMIT, the rest of the code runs.
-	NN is the number of nodes for the smallest arclength. The number of panels is decided based on this number. This code runs a wrapper: NN_start is run and then NN_start+1 is run. If that results in a change in the coefficients of less than tolerance (defined in MOEA_control), that is reported as the answer. Otherwise, NN_start+2 is compared to NN_start+1 (etc) until convergence is achieved.
o	To determine if there are too many panels for WAMIT to run, ‘check_num_panels’ is run
-	‘general_WAMIT_control’ is used to actually run WAMIT

dominance_btwn_two_individuals.py
-	This check whether one individual dominates another

check_num_panels.py
-	This checks how many panels will be needed for WAMIT, given a certain NN and vector. See ‘make_gdf_w_top’ for more details on the mesh

general_WAMIT_control.py
-	This is the main control file for running WAMIT and getting appropriate hydrodynamic coefficients out of it 
-	It calls ‘calculate_volume’ to find mass (assumes freely floating)
-	It also calls ‘delete_WAMIT_files’ which deletes old WAMIT files (because otherwise an error would come up)
-	It calls ‘make_gdf_w_top’ to make the ‘.gdf’ file, ‘make_pot’ to make the ‘.pot’ file, and ‘make_frc_altform_1’ or ‘make_frc’ to make the ‘.frc’ file

delete_WAMIT_files.py
-	This deletes old WAMIT files to avoid an error

make_gdf_w_top.py
-	This code makes the ‘.gdf’ file (geometry file for WAMIT) from the vector, turning the basis functions into a mesh
-	It goes through each of the 5 faces of the WEC and makes the mesh based on the number of nodes on each face (which is determined by the smallest arclength and such that the panels are close to squares)

make_pot.py
-	This makes the ‘.pdf’ file for WAMIT, based on the frequencies and directions of interest and the water depth and location of the fixed point

make_frc.py
-	This makes the ‘.frc’ file for WAMIT when IALTFORM=2

make_frc_alform_1.py
-	code makes the ‘.frc’ file for WAMIT, which specified certain options for outputs etc. 
-	Note that we just used WAMIT to find the hydrodynamic coefficients, not the motions etc 
