%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [T] = solve_fick_problem_1D(x, t, T0, Tw, MatParams)
%SOLVE_FICK_PROBLEM_1D Solves a semi-infinite Fick-problem at a point x at 
%time t with wall temperature Tw  and initial temperature T0 for a 
%given material with parameters MatParams

%% Solve

% get alpha
alpha = MatParams.k / ( MatParams.rho * MatParams.cp );

% get eta
eta = x / ( 2 * sqrt(alpha*t) );

% compute erf(eta)
erfEta = erf(eta);

% compute result  
T = Tw + (T0 - Tw) * erfEta;


end


