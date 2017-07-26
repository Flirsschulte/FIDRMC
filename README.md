# FIDRMC
    Fuzzy Impulse Noise Detection and Reduction Method for Colour Images


    Schulte Stefan, De Witte Valérie, Nachtegael Mike, Van Der Weken Dietrich and Kerre Etienne E. A Fuzzy Two-step Filter for Impulse Noise Reduction From Colour Images. IEEE Transactions on Image Processing 15 (11) (2006) 3567 - 3578 

    Download the zip file here.


    Unzip the file into a certain directory.


    You have to add the folder to the current MATLAB search path by using the commando: addpath('path'), with ‘path’ the directory followed by the new map '/FIDRMC' which should be created.


    The FIDRMC filter is executed by using the Matlab commando:

    Out = FIDRMC(Noise,Org);
    where Noise is the noisy input colour image and Org is the noise-free colour original image.


    In order to overcome some problems of Matlab you have to know that we have implemented most files in C-code. The switching between Matlab and C in realized by the 'mex' (see Matlab help for more information). If you never have used mex files before you probably have to do the following thinks:

    cd 'C:\Directory\'
    the current directory has to be changed to that map that was created in the first step

    mex –setup
    => Choose a compiler

    mex Detection.c
    mex Detectionb.c
    mex DenoiseFIDRMC.c
    mex DenoiseFIDRMC2.c

