%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ Indices ] = convert_from_condensed_index( aMultiIndex, aNumIndices, aLimit )
%CONVERT_FROM_CONDENSED_INDEX Convert condensed multi-index back to individual indices
%   aMultiIndex - multi-index to be converted
%   aNumIndices - number of individual indices to be generated
%   aLimit      - maximum value any individual index may take
%
% OUT:
%   Indices - array of individual indices

%% Process and Check Inputs

% get the max multi-index
MaxMultiIndex = aLimit^aNumIndices;

% check  this value is not exceeded
assert( aMultiIndex <= MaxMultiIndex && aMultiIndex > 0, ...
    'convert_from_condensed_index() - aMultiIndex provided is out of bounds.' );

% initialize array of individual indices
Indices = zeros( 1, aNumIndices );

% convert the multi-index to be 0-based
aMultiIndex = aMultiIndex - 1;

%% Convert

for iInd = aNumIndices:-1:1

    Remainder = mod( aMultiIndex, aLimit ^ (iInd-1) );
    Indices( iInd ) = ( aMultiIndex - Remainder ) / ( aLimit ^ (iInd-1) );
    aMultiIndex = Remainder;

end

% reverse the order of the individual indices
Indices = flip( Indices );

% convert indices to be 1-based
Indices = Indices + 1;

%% return Indices
end

