# Persistent-Homology-Tomato-Leaf-Root
This is the code to compute persistent homology to measure tomato leaf shape, leaf serrations, and root architecture.
*************************************************************************************************************
This package contains all Matlab functions needed to calculate the persistence barcode, bottleneck distance for leaf shape, Euler characteristic curve for serrations, and Betti 1 curve for root architecture used in the analysis of quantitative trait locus (QTL) mapping described in:
Mao Li, Margaret H. Frank, Viktoriya Coneva, Washington Mio, Daniel H. Chitwood, Christopher N. Topp (2018) The persistent homology mathematical framework provides enhanced genotype-to-phenotype associations for plant morphology, Plant Physiology, doi: https://doi.org/10.1104/pp.18.00104

http://www.plantphysiol.org/content/early/2018/06/04/pp.18.00104

Code is implemented by Mao Li (maoli0923@gmail.com).

A few remarks:

(a) The study was based on binary images.

(b) Computing persistence barcode needs to install Javaplex http://appliedtopology.github.io/javaplex/

    (A reduced computational version which does not require to install Javaplex will be released in another following paper.)

(c) To simplify use, create two folders, say, 'Leaves' and 'Roots'. Place all contents of this package and all leaf images in folder 'Leaves' and place all contents of this package and all root images in folder 'roots'.

(d) Images can be in any standard format such as .jpg, .png, .tiff, etc. The default format for leaf image is .jpg, for root image is .png.

(d) The package contains three main scripts (.m files) that perform the following operations: (i) compute persistent barcodes and pairwise bottleneck distance to measure leaf shape; (ii) compute Euler characteristic curve to measure serrations; (iii) calculate the Betti 1 curve to quantify complexity of root architecture.

*************************************************************************************************************************************
First, in Matlab, change Matlab's "Current folder" to the directory 'matlab_example' under 'javaplex' and excute the Matlab command:

    >> load_javaplex

_________________________________________________________________________________________________________________
###1. Measure leaf shape

Then, set the working directory to the folder 'Leaves', and execute the following Matlab command:
 
    >> PersistentHomology_leaf


The output BD_leaf.mat is N-by-N matrix, where N is the number of images of leaves, and BD_leaf(i,j) is the distance between leaf i (with name file(i).name) and leaf j (with name file(j).name).


________________________________________________________________________________________________________________
###2. Measure serrations

Still set the working directory to the folder 'Leaves', and execute the following Matlab command:
 
    >> PersistentHomology_serration

The output serration.mat is N-by-M matrix, where N is the number of images of leaves, and each row is the M-dimensional features to describe the serrations. In this study, the default M=112.

_____________________________________________________________________________________________________________________
###3. Measure 2D projection of root architecture

Finally, set the working directory to the folder 'Roots', and execute the following Matlab command:
     
      >> PersistentHomology_root

The output H1.mat is N-by-M matrix, where N is the number of images of roots, and each row is the M-dimensional features to describe complexity of root branching structure. In this study, the default M=80.
