%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [T] = solve_stefans_problem_1D(x, t, Tw, MatParams)
%SOLVE_STEFANS_PROBLEM_1D Solves for the temperatur the 1D phase change 
%Stefan's problem at a point x at time t with wall temperature Tw for a 
%given material with parameters MatParams

%% Solve

% get alpha
alpha = MatParams.k / ( MatParams.rho * MatParams.cp );

% get eta
eta = x / ( 2 * sqrt(alpha*t) );

% compute erf(eta)
erfEta = erf(eta);

% compute Stefan's number
Ste = ( MatParams.cp / MatParams.Lh ) * (Tw - MatParams.Tm);

% solve for beta
beta_eqn = @(beta) ( beta * exp(beta^2) * erf(beta) - 0.5 * Ste );
beta = fzero( beta_eqn, 1.0 );

% compute erf(beta)
erfBeta = erf(beta);

% get position of melting front
delta = 2 * eta * sqrt(alpha * t);

% debug
% fprintf( 'Melting Front at: %f \n x = %f, t = %f \n', delta, x, t);

% compute result  
T = Tw + (MatParams.Tm - Tw) * (erfEta/erfBeta);

% clean results behind melting front
if (T < MatParams.Tm)
    T = MatParams.Tm;
end  

end


