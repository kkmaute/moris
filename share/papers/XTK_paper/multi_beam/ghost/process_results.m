%% Script to check rank of the linear system
close all;
clear;
clc;

%% Input and Output values

% offset values the simulations were performed with
OffSetValsStr = ["-140" "-139" "-120" "-101" "-100" "-080" "-060" "-040" "-039" "-020" "-001" "000" "020" "040" "060" "061" "080" "099" "100" "120" "140" "160" "161" "180" "199" "200"];
NumVals = length( OffSetValsStr );

% initialize output values
OffSetVals = zeros( NumVals, 1 );
CondNumsGhostOff = zeros( NumVals, 1 );
CondNumsGhostOn = zeros( NumVals, 1 );
NumDofsGhostOff = zeros( NumVals, 1 );
NumDofsGhostOn = zeros( NumVals, 1 );

%% Loop reading and processing results for all simulations

for i = 1:NumVals
    
    % get the offset value
    OffSetStr = OffSetValsStr( i );
    OffSetVal = str2double( OffSetStr ) / 10.0;
    OffSetVals( i ) = OffSetVal;

    % get the condition number for Ghost OFF
    Dir = "Ghost_Off/OffSet_" + OffSetStr;
    cd( Dir );
    load( 'SOE.0.jac.dat' );
    SparseJacobian = sparse( SOE_0_jac(:,1), SOE_0_jac(:,2), SOE_0_jac(:,3) );
    Jacobian = full( SparseJacobian );
    CondNumsGhostOff( i ) = cond( Jacobian );
    NumDofsGhostOff( i ) = length( Jacobian );
    cd ../../

    % get the condition number for Ghost ON
    Dir = "Ghost_On/OffSet_" + OffSetStr;
    cd( Dir );
    load( 'SOE.0.jac.dat' );
    SparseJacobian = sparse( SOE_0_jac(:,1), SOE_0_jac(:,2), SOE_0_jac(:,3) );
    Jacobian = full( SparseJacobian );
    CondNumsGhostOn( i ) = cond( Jacobian );
    NumDofsGhostOn( i ) = length( Jacobian );
    cd ../../

end

%% Plotting

figure;
semilogy( OffSetVals, CondNumsGhostOff, "-square", 'LineWidth', 2, 'MarkerSize', 5 );
hold on;
semilogy( OffSetVals, CondNumsGhostOn, "-o", 'LineWidth', 2, 'MarkerSize', 5 );
xlabel( "o (x-Offset)" );
ylabel( "cond(A)" );
legend( 'Ghost: Off', 'Ghost: On', "Location", "northwest" );

figure;
plot( OffSetVals, NumDofsGhostOn, "-square", 'LineWidth', 2, 'MarkerSize', 5 );
hold on;
plot( OffSetVals, NumDofsGhostOff, "-o", 'LineWidth', 2, 'MarkerSize', 5 );
xlabel( "o (x-Offset)" );
ylabel( "n_{DOF}" );
legend( 'Ghost: Off', 'Ghost: On', "Location", "northwest" );




