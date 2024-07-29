%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Semi-Analytical Solution for a 1D Stefan's Problem
close all;
clear;
clc;

%% Material Parameters

% thermal conductivity
MatParams.k = 0.21;

% density
MatParams.rho = 750;

% heat capacity
MatParams.cp = 2.4e3;

% latent heat
MatParams.Lh = 175.0e3;

% melting temperature
MatParams.Tm = 313.0;

%% Problem Setup (Geometry & BCs)

% wall temperature
Tw = 350.0;

% length of interest
lengthBar = 0.12;

% number of elements for evaluation
nElements = 100;

% time frame of solution
tmax = 16 * 3600;

% number of timesteps for evaluation
nTimeSteps = 1;

%% Compute Solution

% initialize
deltat = tmax / nTimeSteps;
deltax = lengthBar / nElements;
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
        SolutionT(it,iNode) = solve_stefans_problem_1D(x, t, Tw, MatParams);
               
    end
    
end

%% Plot Solution

% read numerical solution
% [xNodesNumericalGGLS, NumericalTGGLS] = read_solution_bar('Comsol_16h00.csv', 4, 1);
% [xNodesNumerical, NumericalT] = read_solution_bar('Comsol_wo_GGLS_16h00.csv', 4, 1);

% initialize vector of node points
xNodes = linspace(0,lengthBar,nElements+1);

% for each time step
for it = nTimeSteps+1:nTimeSteps+1
    
    % get point in time
    t = (it - 1) * deltat;
    
    % create figure and label
    figure;
    hold on;
    title(['Temperature solution at t = ',num2str(t/3600),' hours.'])
    axis([0, lengthBar, MatParams.Tm , Tw]);
    
    % plot solution
    plot( xNodes, SolutionT(it,:), 'k', 'LineWidth', 1.6 );
    
    % plot numerical solution
%     plot( xNodesNumerical, NumericalT, 'b', 'LineWidth', 1.6 );
%     plot( xNodesNumericalGGLS, NumericalTGGLS, 'r', 'LineWidth', 1.6);
    
    % add legend
%     legend('analytic','numeric', 'numeric /w GGLS');
    ylabel('Temperature');
    xlabel('x');
           
end

