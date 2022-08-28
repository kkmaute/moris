%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Read matrices

load matrices_rho.dat
matrices_rho = full( spconvert( matrices_rho ) );
load matrices_vel.dat
matrices_vel = full( spconvert( matrices_vel ) );
load matrices_tmp.dat
matrices_tmp = full( spconvert( matrices_tmp ) );
load matrices.dat
matrices = full( spconvert( matrices ) );
load matrices_left_nitsche.dat
matrices_left_nitsche = full( spconvert( matrices_left_nitsche ) );


%% Compare matrices

figure;
spy(matrices_rho);
title('only rho contribution');
condA_rho = cond(matrices_rho)

figure;
spy(matrices_vel);
title('only vel contribution');
condA_vel = cond(matrices_vel)

figure;
spy(matrices_tmp);
title('only tmp contribution');
condA_tmp = cond(matrices_tmp)

