# Design of Experiments on Non-Convex Regions
Code and dataset for our 2017 Technometrics [paper](https://drive.google.com/file/d/1UKfKGxTT_hBuesyVI3cJORo3OEzXwMxa/view?usp=drive_link) "Design and Analysis of Experiments on Non-Convex Regions". Each folder contains its own "ReadMe.pdf" file with explicit instrcuctions on running its own designated code.

## 2D Horseshoe
This folder contains the R code and the datasetsused for the 2D Horseshoe example appearing in

Sections 3 and 4 of the paper. If you wish to reproduce the simulation results, we recommend that

you run the different code parts in the following order:

1. **Horseshoe_Pre_Optimization.R:** This file sources the functions file *Pre_Optimization_Functions.R*, receives grid density (epsilon) and embedding dimension (numapproxevals) as arguments and writes the following files into the /Simulation/Auxiliary_Files directory:

 - *Horseshoe_Candidates_(epsilon).txt:* list of candidates on a grid covering the horseshoe in the original Euclidean coordinates.
 - *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt:* the original candidates in their representation in the new coordinate system in the embedding space after MDS.
 - *Horseshoe_Geodesic_Distances_(epsilon).txt:* The approximate geodesic distance matrix for thecandidates.
 - *Horseshoe_Embedded_Candidates_Full_(epsilon).txt:* the original candidates in their representation in the new coordinate system in the embedding space after MDS, including ALL the eigenvectors corresponding to non­‐negative eigenvalues. Serves to find the best embedding dimension.
   
2. **Horseshoe_Maximin_Design.R:** This file reads the geodesic distance matrix *Horseshoe_Geodesic_Distances_(epsilon).txt* as well as the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Maximin geodesic distance criterion. The results are written into the */Simulation/Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

3. ­**Horseshoe_Euclidean_Maximin_Design.R:** this file reads the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size design.size using an exchange algorithm, based on the Maximin Euclidean distance criterion. The results are written into the */Simulation/Euclidean_Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_Euclidean_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_Euclidean_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

4. ­**Horseshoe_IMSPE_Designs.R:** The code sources *Pre_Optimization_Functions.R*, the simulation design file */Simulation/Simulation_Design_List.txt*, the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Minimum IMSPE criterion. The results are written into the /Simulation/IMSPE_Designs directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_IMSPE_Design_(#).txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_IMSPE_Design_(#).txt:* The resultant design in the new coordinate system in the embedding space.

5. ­IMSPE_Comparison.R : this part of the code reads the output of the previous ones, and produces the results of the simulation study of Section 4. Along the way different tables are printed to the screen: **Maximin.Eff** contains a detailed account of the relative efficiency of each of the maximin (geodesic and Euclidean distance) design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs. Likewise, **IMSPE.all** contains a detailed account of the relative efficiency of each of the IMSPE­‐optimal design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs.

## 3D Horseshoe
This folder contains the R code and the datasetsused for the 2D Horseshoe example appearing in

Sections 3 and 4 of the paper. If you wish to reproduce the simulation results, we recommend that

you run the different code parts in the following order:

1. **Horseshoe__3D_Pre_Optimization.R:** This file sources the functions file *Pre_Optimization_Functions.R*, receives grid density (epsilon) and embedding dimension (numapproxevals) as arguments and writes the following files into the /Simulation/Auxiliary_Files directory:

 - *Horseshoe_3D_Candidates_(epsilon).txt:* list of candidates on a grid covering the horseshoe in the original Euclidean coordinates.
 - *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt:* the original candidates in their representation in the new coordinate system in the embedding space after MDS.
 - *Horseshoe_3D_Geodesic_Distances_(epsilon).txt:* The approximate geodesic distance matrix for thecandidates.
 - *Horseshoe_3D_Embedded_Candidates_Full_(epsilon).txt:* the original candidates in their representation in the new coordinate system in the embedding space after MDS, including ALL the eigenvectors corresponding to non­‐negative eigenvalues. Serves to find the best embedding dimension.
   
2. **Horseshoe_3D_Maximin_Design.R:** This file reads the geodesic distance matrix *Horseshoe_3D_Geodesic_Distances_(epsilon).txt* as well as the set of candidates *Horseshoe_3D_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Maximin geodesic distance criterion. The results are written into the */Simulation/Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_3D_Size_(design.size)_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

- Later, an optimal subset of size design.size is chosen using an exchange algorithm, based on the Maximin Euclidean distance criterion. The results are written into the */Simulation/Euclidean_Maximin_Designs* directory, in the form of the following files: 

 - *Horseshoe_3D_Size_(design.size)_Euclidean_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_Euclidean_Maximin_Design.txt:* the resultant design in the new coordinate system in the embedding space.

3. ­**Horseshoe_IMSPE_Designs.R:** The code sources *Pre_Optimization_Functions.R*, the simulation design file */Simulation/Simulation_Design_List.txt*, the set of candidates *Horseshoe_3D_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Minimum IMSPE criterion. The results are written into the /Simulation/IMSPE_Designs directory, in the form of the following files:

 - *Horseshoe_3D_Size_(design.size)_IMSPE_Design_(#).txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_IMSPE_Design_(#).txt:* The resultant design in the new coordinate system in the embedding space.

. ­IMSPE_Comparison.R : this part of the code reads the output of the previous ones, and produces the results of the simulation study of Section 4. Along the way different tables are printed to the screen: **Maximin.Eff** contains a detailed account of the relative efficiency of each of the maximin (geodesic and Euclidean distance) design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs. Likewise, **IMSPE.all** contains a detailed account of the relative efficiency of each of the IMSPE­‐optimal design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs.
