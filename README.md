# PresolarSiCclustering2021

This R code is meant to be used for clustering data on presolar silicon carbide (SiC) grains. DOI: 10.5281/zenodo.4304837

It is linked to the article currently in revision for publication in The Astrophysical Journal Letters:
"Cluster analysis of presolar silicon carbide grains: evaluation of their classification and astrophysical implications"
by Asmaa Boujibar1, Samantha Howell1,2, Shuang Zhang1, Grethe Hystad3, Anirudh Prabhu4, Nan Liu5, Thomas Stephan6,7, Shweta Narkar4, Ahmed Eleish4, Shaunna M. Morrison1, Robert M. Hazen1, Larry R. Nittler1.

Affiliations: 1Earth and Planets Laboratory, Carnegie Institution for Science, 2Washington College, 
3Purdue University Northwest, 4Tetherless World Constellation, Rensselaer Polytechnic Institute, 
5Washington University in St. Louis, Sciences, 6The University of Chicago, 7Chicago Center for Cosmochemistry


Original data (excel file PGD_SiC_2020-08-18.xlsx) is open access and can be found on the website of Washington University in St. Louis:
https://presolar.physics.wustl.edu/presolar-grain-database/

In this code, the chosen clustering technique is a model-based clustering algorithm assuming mixtures of Gaussian distributions: Mclust R package (Scrucca et al. 2016).
We chose in this example the four features: 12C/13C, 14N/15N, δ29Si/28Si and δ30Si/28Si ratios. 

An excel file is linked to this code, including grain data and clustering results.
Discussions of the results and their astrophysical implications will be available in the article once published by the journal The Astrophysical Journal Letters.
