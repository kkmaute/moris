%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%



%% #######################   Header Start ################################
%
% File Name:            CompareFDSweepHDF5.m
% Author:               Sam Gates
% Date Created:         October 12, 2021
%
%*************************************************************************
%
% Description:   
%      This script is meant to open a given HDF5 file of a FDsweep and to
%      compare the analytical solution to the central, foward, and backward
%      FD calculations.
% 
% Inputs:
%      - tFileName:     File name pattern to open without the file type,
%                       eg, "test_20.1.hdf5"
%
% Outputs:
%      - Max difference
%      - Max difference plot
%
%  ######################   Header Stop ##################################

%% ######################   Code Start  ##################################
clear;
clc;
close all;

tFileName = "test.hdf5";  %Filename
    
% Open the file
tFileID = H5F.open(tFileName,'H5F_ACC_RDONLY','H5P_DEFAULT');
i = 1;

% open data, send to array, close for all 8 data sets
tDataSetID = H5D.open(tFileID,'constraint_gradients eval_1-1 analytical');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'constraint_gradients eval_1-1 epsilon_1-1 fd_backward');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'constraint_gradients eval_1-1 epsilon_1-1 fd_central');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'constraint_gradients eval_1-1 epsilon_1-1 fd_forward');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'objective_gradients eval_1-1 analytical');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'objective_gradients eval_1-1 epsilon_1-1 fd_backward');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'objective_gradients eval_1-1 epsilon_1-1 fd_central');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

tDataSetID = H5D.open(tFileID,'objective_gradients eval_1-1 epsilon_1-1 fd_forward');
tData(:,i) = H5D.read(tDataSetID);
H5D.close(tDataSetID);
i = i + 1;

% close hdf5 file
H5F.close(tFileID);

% ######################   Code Stop  ##################################

%% #####################   Plotting Start  #############################

% plotting constraint error
figure(1)
tiledlayout(3,1)

nexttile
plot( abs( tData(:,1) - tData(:,2) ) )
title('Constraint Gradients FD Backward')

nexttile
plot( abs( tData(:,1) - tData(:,3) ) )
title('Constraint Gradients FD Central')

nexttile
plot( abs( tData(:,1) - tData(:,4) ) )
title('Constraint Gradients FD Forward')

saveas(figure(1),'Constraint_Difference','png')

% plotting objective error
figure(2)
tiledlayout(3,1)

nexttile
plot( abs( tData(:,5) - tData(:,6) ) )
title('Objective Gradients FD Backward')

nexttile
plot( abs( tData(:,5) - tData(:,7) ) )
title('Objective Gradients FD Central')

nexttile
plot( abs( tData(:,5) - tData(:,8) ) )
title('Objective Gradients FD Forward')

saveas(figure(2),'Objective_Difference','png')

% ######################   Plotting Stop  ################################





