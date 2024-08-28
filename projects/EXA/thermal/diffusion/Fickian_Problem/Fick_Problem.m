%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Semi-Analytical Solution for a 1D Stefan's Problem
%close all;
clear;
clc;

%% Material Parameters

% thermal conductivity
MatParams.k = 2.1e-7;

% density
MatParams.rho = 0.750;

% heat capacity
MatParams.cp = 2.4e0;

%% Problem Setup (Geometry & BCs)

% initial temperature
T0 = 313.0;

% wall temperature
Tw = 350.0;

% length of interest
length = 0.25;

% number of elements for evaluation
nElements = 50;

% time frame of solution
tmax = 16.0/3.0 * 3600.0;

% number of timesteps for evaluation
nTimeSteps = 1;

%% Compute Solution

% initialize
deltat = tmax / nTimeSteps;
deltax = length / nElements;
SolutionT = zeros(nTimeSteps+1,nElements+1);

% for each time step
for it = 1:nTimeSteps+1
    
    % get point in time
    t = (it - 1) * deltat;
    
    % for each node
    for iNode = 1:nElements+1
        
        % get point in space
        x = (iNode - 1) * deltax;
        
        % compute solution at x & t
        SolutionT(it,iNode) = solve_fick_problem_1D(x, t, T0, Tw, MatParams);
               
    end
    
end

%% Plot Solution

% read numerical solution
% [xNodesNumerical, NumericalT] = read_solution_bar('Fick_5h20.csv', 4, 1);

% initialize vector of node points
xNodes = linspace(0,length,nElements+1);

% for each time step
for it = nTimeSteps+1:nTimeSteps+1
    
    % get point in time
    t = (it - 1) * deltat;
    
    % create figure and label
    figure;
    hold on;
    title(['Temperature solution at t = ',num2str(t/3600),' hours.'])
    axis([0, length, T0 , Tw]);
    
    % plot solution
    plot( xNodes, SolutionT(it,:), '-x' );
    
    % plot numerical solution
%     plot( xNodesNumerical, NumericalT );
    
    % add legend
%     legend('analytic','numeric');
           
end

