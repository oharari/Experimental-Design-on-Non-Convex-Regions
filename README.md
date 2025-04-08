# Design of Experiments on Non-Convex Regions
Code and dataset for our [2017 Technometrics paper](https://www.tandfonline.com/doi/10.1080/00401706.2015.1115674) (also available as a [preprint](https://drive.google.com/file/d/1UKfKGxTT_hBuesyVI3cJORo3OEzXwMxa/view?usp=drive_link)) "Design and Analysis of Experiments on Non-Convex Regions". Each folder contains its own "ReadMe.pdf" file with explicit instrcuctions on running its own designated code.

## 2D Horseshoe
This folder contains the R code and the datasetsused for the 2D Horseshoe example appearing in Sections 3 and 4 of the paper. If you wish to reproduce the simulation results, we recommend that you run the different code parts in the following order:

1. **Horseshoe_Pre_Optimization.R:** This file sources the functions file *Pre_Optimization_Functions.R*, receives grid density (epsilon) and embedding dimension (numapproxevals) as arguments and writes the following files into the /Simulation/Auxiliary_Files directory:

 - *Horseshoe_Candidates_(epsilon).txt:* A list of candidates on a grid covering the horseshoe in the original Euclidean coordinates.
 - *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt:* The original candidates in their representation in the new coordinate system in the embedding space after MDS.
 - *Horseshoe_Geodesic_Distances_(epsilon).txt:* The approximate geodesic distance matrix for thecandidates.
 - *Horseshoe_Embedded_Candidates_Full_(epsilon).txt:* The original candidates in their representation in the new coordinate system in the embedding space after MDS, including ALL the eigenvectors corresponding to non­‐negative eigenvalues. Serves to find the best embedding dimension.
   
2. **Horseshoe_Maximin_Design.R:** This file reads the geodesic distance matrix *Horseshoe_Geodesic_Distances_(epsilon).txt* as well as the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Maximin geodesic distance criterion. The results are written into the */Simulation/Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

3. ­**Horseshoe_Euclidean_Maximin_Design.R:** This file reads the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size design.size using an exchange algorithm, based on the Maximin Euclidean distance criterion. The results are written into the */Simulation/Euclidean_Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_Euclidean_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_Euclidean_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

4. ­**Horseshoe_IMSPE_Designs.R:** The code sources *Pre_Optimization_Functions.R*, the simulation design file */Simulation/Simulation_Design_List.txt*, the set of candidates *Horseshoe_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Minimum IMSPE criterion. The results are written into the /Simulation/IMSPE_Designs directory, in the form of the following files:

 - *Horseshoe_Size_(design.size)_IMSPE_Design_(#).txt:* The resultant design in the original coordinates.
 - *Horseshoe_Size_(design.size)_Embedded_IMSPE_Design_(#).txt:* The resultant design in the new coordinate system in the embedding space.

5. ­**IMSPE_Comparison.R:** This part of the code reads the output of the previous ones, and produces the results of the simulation study of Section 4. Along the way different tables are printed to the screen: **Maximin.Eff** contains a detailed account of the relative efficiency of each of the maximin (geodesic and Euclidean distance) design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs. Likewise, **IMSPE.all** contains a detailed account of the relative efficiency of each of the IMSPE­‐optimal design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs.

## 3D Horseshoe
This folder contains the R code and the datasetsused for the 3D Horseshoe example appearing in Sections 4 of the paper. If you wish to reproduce the simulation results, we recommend that you run the different code parts in the following order:

1. **Horseshoe__3D_Pre_Optimization.R:** This file sources the functions file *Pre_Optimization_Functions.R*, receives grid density (epsilon) and embedding dimension (numapproxevals) as arguments and writes the following files into the /Simulation/Auxiliary_Files directory:

 - *Horseshoe_3D_Candidates_(epsilon).txt:* A list of candidates on a grid covering the horseshoe in the original Euclidean coordinates.
 - *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt:* The original candidates in their representation in the new coordinate system in the embedding space after MDS.
 - *Horseshoe_3D_Geodesic_Distances_(epsilon).txt:* The approximate geodesic distance matrix for thecandidates.
 - *Horseshoe_3D_Embedded_Candidates_Full_(epsilon).txt:* The original candidates in their representation in the new coordinate system in the embedding space after MDS, including ALL the eigenvectors corresponding to non­‐negative eigenvalues. Serves to find the best embedding dimension.
   
2. **Horseshoe_3D_Maximin_Design.R:** This file reads the geodesic distance matrix *Horseshoe_3D_Geodesic_Distances_(epsilon).txt* as well as the set of candidates *Horseshoe_3D_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Maximin geodesic distance criterion. The results are written into the */Simulation/Maximin_Designs* directory, in the form of the following files:

 - *Horseshoe_3D_Size_(design.size)_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

- Later, an optimal subset of size design.size is chosen using an exchange algorithm, based on the Maximin Euclidean distance criterion. The results are written into the */Simulation/Euclidean_Maximin_Designs* directory, in the form of the following files: 

 - *Horseshoe_3D_Size_(design.size)_Euclidean_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_Euclidean_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

3. ­**Horseshoe_IMSPE_Designs.R:** The code sources *Pre_Optimization_Functions.R*, the simulation design file */Simulation/Simulation_Design_List.txt*, the set of candidates *Horseshoe_3D_Candidates_(epsilon).txt* and its respective list of embedded points *Horseshoe_3D_Embedded_Candidates_(numapproxevals)Dim_(epsilon).txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Minimum IMSPE criterion. The results are written into the /Simulation/IMSPE_Designs directory, in the form of the following files:

 - *Horseshoe_3D_Size_(design.size)_IMSPE_Design_(#).txt:* The resultant design in the original coordinates.
 - *Horseshoe_3D_Size_(design.size)_Embedded_IMSPE_Design_(#).txt:* The resultant design in the new coordinate system in the embedding space.

4. ­**IMSPE_Comparison.R:** This part of the code reads the output of the previous ones, and produces the results of the simulation study of Section 4. Along the way different tables are printed to the screen: **Maximin.Eff** contains a detailed account of the relative efficiency of each of the maximin (geodesic and Euclidean distance) design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs. Likewise, **IMSPE.all** contains a detailed account of the relative efficiency of each of the IMSPE­‐optimal design with respect to the true parameters of the Gaussian process model, while **Ave.Maximin.Eff** contains the average relative efficiency of the two for each number of runs.

# Glacier

This folder contains the R code and the datasets used for theglacier application appearing in Section 5 of the paper. If you wish to reproduce the results, we recommend that you run the different code steps in the order they appear. However, you should be warned that the entire process might require over a day, therefore we suggest that you save the different folders in the proposed file architecture, and run only parts 3 and 7 (maximin geodesic distance and IMSPE­‐optimal designs, respectively). Hereby is a short description of the different steps:

1. **Glacier_Scaling.R:** Reads the original locations on the glacier from the file *Glacier_Interior_and_Boundary.txt* (in the */Auxiliary_Files directory*, first 247 rows define the glacier’s boundary), scales and centers them so that they are contained in the unit cube and writes the new locations into the *Glacier_Interior_and_Boundary_Scaled.txt* file in the same directory.

2. **Pre_Optimization.R:** sources the functions file *Pre_Optimization_Functions.R*, performs Isomap on the glacier data (using a 2X10^­‐5 ball about each location to define its neighbors) and writes the following files into the */Auxiliary_Files* directory:

 - *Glacier_Embeddings_4D.txt:* The original locations in their representation in the new coordinate system in the embedding space after MDS.
 - *Glacier_Geodesic_Distances.txt:* The approximate geodesic distance matrix for the candidates.

3. **Glacier_Maximin_Design.R:** Reads the set of scaled locations *Glacier_Interior_and_Boundary_Scaled.txt* and its respective list of embedded locations *Glacier_Embeddings_4D.txt* and chooses an optimal subset of size *design.size* using an exchange algorithm, based on the Maximin geodesic distance criterion. The results are written into the */Glacier_Maximin_Designs* directory, in the form of the following files:

 - *Glacier_Size_(design.size)_Maximin_Design.txt:* The resultant design in the original coordinates.
 - *Glacier_Size_(design.size)_Embedded_Maximin_Design.txt:* The resultant design in the new coordinate system in the embedding space.

4. **Embedding_Integration_Pts.R:** Reads a set of 5922 locations on the glacier from the file *Glacier_Integration_Points.txt* in the */Auxiliary_Files* directory and uses the Nystrom approximation to embed those in the higher dimensional space, writing the output into *Embedded_Integration_Pts_Glacier_Scaled_all.txt* in the same directory.

5. **Glacier_Uniform_Integration_Pts.R:** Reads *Embedded_Integration_Pts_Glacier_Scaled_all.txt* (written in the previous step) and selects a subset of uniformly spaced location for integration over the embedded region, using the maximin geodesic distance (more precisely, its approximation) criterion. The selected subset is written into the file *Glacier_Uniform_Integration_Pts.txt* in the */Auxiliary_Files* directory.

6. ­**Embedding_Glacier_Data_2009.R:** Reads the 2009 glacier ablation data set (consisting of 14 observations) from the *Glacier_Data_2009.txt* file in the */Glacier_2009_Measurements* directory, scales it so that it is contained in the unit cube (in a way that is consistent with the scaling of the candidates points), uses Nystrom’s method to embed it in the 4­‐dimensional region and writes the embedded data set into the *Embedded_2009_Data_Glacier_Scaled.txt* file in the same directory.

7. **Glacier_IMSPE_Optimal_Design.R:** Sources the functions file *IMSPE_Optimal_Design_Functions.R* and generates an IMSPE­‐optimal design of any size on the glacier, based on the geodesic distance. Uses the maximin geodesic distance design of the same size, generated in step 3 and the maximum likelihood estimators of the Gaussian process parameters, estimated from the 2009 data. The results are written into *Glacier_Size_(design.size)_Embedded_IMSPE_Design_2009.txt* and *Glacier_Size_(design.size)_IMSPE_Design_2009.txt* in the */Glacier_IMSPE_Optimal_Designs/2009_Designs* directory.
