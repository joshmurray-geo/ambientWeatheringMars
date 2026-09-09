Repo contains the following:
dryPerplexInput.dat – an example Perple_X input file used to generate our subsequent compositions and mineral assemblages. 

dryPerplex_minModes.tab, solidSolution.tab, and totalComp.tab are those outputs. Mineral abundances, solid solution 
composition, and bulk composition of the system, respectively. 

for the rate calculation, we include the following data: thompsonLoringWaterFilms.xlsx, oconnor04appendix.xls, tempParameters.mat,
and mlParameters.mat. The first two are data from Thompson et al., 2024 and O'Connor et al., 2004 (see text). The .mat files are 
vectors of fit parameters used to extrapolate for our rate calculation. 

finally, we include various MATLAB scripts used for calculation and synthesis of the data. plumx.m and labelPolygon.m are used for
graphing pseudosections in MATLAB based upon mineral modes. examplePseudosect.m is a script that calls those graphing functions in
order to produce the dry pseudosection. And temperatureAnalysis.m and monolayerAnalysis.m yield the aforementioned fit parameters 
and uncertainties. 
