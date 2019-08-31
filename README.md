# PhyloDon
    A series of R functions that allow users to "Roll back" phylogenies to earlier points of taxonomic knowledge and track a series of metrics, inclidind Pybus and Harvey's Gamma, Tree length, phylogenetic distance (PD) and average new branch lengths, as well as functions that provide frames of tree-growth through time, and model Pybus and Harvey's gamma forward from previous points in taxonomic knowledge.

# NewTaxNo
    report the number of newly described taxa per year, or block of years, within a phylogeny

# GammaPlot
    plot the gamma values of rolled back tree. Also includes a simulated gamma for each time slice using the tree of same OTU number with a random set of taxa removed.

# LTTTTVideo 
    (Lineages Through Time Through Time) = create a series of sequential images of LTT plots for a range of years

# TreeGrow
reate a series of sequential images of trees for a range of years

# PDPlot 
    plot the Phylogenetic distance of taxa described within a specific time-window, PD = total length of distance matrix: PDcorr= PD divided by number of taxa included in the tree: PDcorr2= PD divided by the total pairwise distance of the total tree.

# TLPlot
    plot the Treelength values of rolled back tree. Also includes a simulated treelength for each time slice using the tree of same OTU numebr with a random set of taxa removed.

# BLplot 
    plot the average and max min lengths of new Branchlengths between sucessive time slices)

# GammaPredict 
    predict the range of gama values of "final tree" run forward from various starting points throughout the taxonomic history 

# GammaRandomPlot 
    this adds X random tips to a tree n, where X is the differece in taxa between n and n+1

# GammaTaxRandomPlot
    This function randomizes the taxonomic positions, and runs multiple Gammaplot itterations, comparing that to the itterated randomly removed specimens.

# GammaModify
    Modify branchlength distribution to match a particular value of Gamma

# lambdaPredict
    Idenitifies the phylogenetic signal (lambda) of X characters, simulated under Brownian motion, across a tree's taxonomic history. 

# RatePredict
    Idenitifies the evolutioniary rate (mean PIC^2) of X characters, simulated under Brownian motion, across a tree's taxonomic history.
  
