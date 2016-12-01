===========================
MATLAB Code for "Solving and Simulating an Unbalanced Growth Model using Linearization about the Current State"
===========================

This repository contains the computational code used in the paper, "Solving and Simulating an Unbalanced Growth Model using Linearization about the Current State," by Kerk L. Phillips.  It also contains a PDF version of the October 2016 working paper and a separate technical appendix which gives a derivation of the Euler errors and the details of the equations which define the simulated models.

MATLAB code
===========
The MATLAB code for this paper is the file UnBal.m in the MATLAB folder. All supporting fuctions and m-files are included in the same directory.  JRstat.m runs the balanced growth model with Jaimovich-Rebelo preferences using a stationarized version of the model.  It uses both the SSL and CSL methods.  GHH.m runs the unbalanced growth model with Greenwood, Hercowitz and Huffmans preferences.  It uses only the CSL method.  Changes to model and other parameter values can be made in these files.  There is no need to change anything in the other m file which are called as functions.