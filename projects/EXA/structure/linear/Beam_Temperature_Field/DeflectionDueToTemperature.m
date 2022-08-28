%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

clear
clc
% This script shows the analytical vertical deflection of a cantilever
% beam due to a temperature gradient from top to bottom
% https://nptel.ac.in/content/storage2/courses/105101085/downloads/lec-26.pdf

% Inputs
E = 1;
nu = 0.3;
CTE = 0.01;
tRefTemp = 0;

% Beam dimensions
L = 20;
d = 1;

% Temperature Fiel
T_top = 1;
T_bot = 0;

% Vertical Deflection
delta_y = CTE * L^2 * (T_top - T_bot) / (2*d)
