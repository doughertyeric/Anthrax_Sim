# Anthrax_Sim

This repository consists of the files required to run and analyze an agent-based model of ungulates moving across a heterogeneous landscape. The landscape is populated with locally-infectious zones (LIZs) that are distributed non-randomly, which represent carcasses of previously infected ungulates (either springbok-sized, zebra-sized, or elephant-sized) that succumbed to the disease. The ultimate purpose of the simulation was to calculate probabilities of contact between agents and these LIZs. For more details on the parameterization of prupose of the simulation, refer to "Using movement data to estimate contact rates in a simulated environmentally-transmitted disease system" (doi: https://doi.org/10.1101/261198).

# File 1: AnthraxSim.R

A file consisting of a set of functions intended to simulate the heterogeneous landscape, the unmoving LIZ sites on that landscape, and the live ungulate agents moving across the landscape. The functions contained are described briefly below:
- a) Basic functions: meant to offer outputs used in several other, more specific, functions (within.func, outside.func, anglefun, movement)
- b) Initialization functions: meant to create a landscape with boundaries defined by the user. Each cell is assigned a forage quality value as well. The edge of the landscape can then be adjusted to decrease the forage quality within a certain distance of the 'fence'. Initial LIZ and agent characteristics and locations are then assigned based on the heterogeneous habitat quality. In addition, the behavioral states of the agents are assigned based on the distribution observed in empirical zebra trajectories. Finally, a green.up function is used to increase the forage quality values of the cells on which LIZs were placed to reflect the observed dynamics associated with vegetation growth in the years following carcass deposition
- c) Agent state processes: meant to alter the states and positions of agents throughout the simulation. move.func.simp is used to designate a particular step length and turning angle based on the current behavioral state of the individual. In state 1 (resting), steps are small and undirected, in state 3 (directed movement), steps are longer and concentrated around a mean direction, and in state 2 (foraging), steps are of intermediate length and are directed towards the highest quality cell within the individual's perceptual range. In addition, the state.vector function shifts the agents behavioral state according to probabilities derived from the empirical paths.
- d) Landscape processes: meant to alter the forage quality throughout the course of the simulation. Both specific (based on the actual agents' positions) and general (based on a somewhat random distribution of untracked agents) feeding occur and decrease habitat quality. Vegetation also grows according to a logistic growth function.
- e) Post run: meant to record the outputs of the simulation. The paths (with 20 minute temporal resolution) are saved for each of the tracked individuals.
The final portion of the code actually runs the simulation over a single anthrax season in the form a series of loops.

# File 2: kmeans_Clustering.R
A file that attempts to place each step (of each path) into one of three behavioral states based on the step size distributions associated with each. It draws upont the kmeans clustering method, and the result is the assignment of each of the points into the foraging state (or one of the other two) for subsequent analyses at broad and fine scales.

# File 3: Broad_Scale_Analysis.R 
A file that uses two alternative methods applied at the broad (home range) scale to estimate probability of contact between agents and LIZs. The first method draws on the minimum convex polygon home range and the second draws upon the LoCoH method (which is optimized based on a modified version of the optimization algorithm implemented in the Movement Ecology article "A cross-validation-based approach for delimiting reliable home range estimates" (Dougherty et al. 2017). The full paths and the reduced paths with only foraging points (based on the kmeans clustering procedure) are analyzed using these approaches.

# File 4: AnthraxSim_STP_Sim.R
A file enabling the simulation of steps in between recorded fixes of each animal. Due to the temporal resolution of the empirical data upon which the ABM is based, the positions of agents are recorded once every 20 minutes. This code draws upon the space-time prism (STP) method to simulate the potential path taken by the animal in between these fixes. By simulating these potential paths numerous times, a probability surface emerges, which forms the basis of the Fine-scale Analysis (File 5)

# File 5: Fine_Scale_Analysis.R 
A file that takes the output from the AnthraxSim_STP_Sim.R file and creates a density raster in order to calculate the probability of an agent being in the same place as a LIZ (based on a generalized LIZ raster layer). The full paths and the reduced paths with only foraging points (based on the kmeans clustering procedure) are analyzed using this approach.

