%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ MultiIndex, MaxMultiIndex ] = convert_to_condensed_index( aIndices, aLimit )
%CONVERT_TO_DOUBLE_INDEX convert individual indices into a condensed multi-index
%   aIndices - array of individual indices to be converted into a multi-index
%   aLimit   - maximum value any individual index may take
%
% OUT:
%   MultiIndex    - multi-index
%   MaxMultiIndex - maximum multi-index

%% Process and Check Inputs

% get the number of indices
NumIndices = length( aIndices );

% check that no index is out of bounds
for iInd = 1:NumIndices
    assert( aIndices( iInd ) <= aLimit &&  aIndices( iInd ) > 0, ...
        'convert_to_condensed_index() - index #%i exceeds bounds: %i > %i', ...
        iInd, aIndices( iInd ), aLimit );
end


%% Convert

% get the max multi-index
MaxMultiIndex = aLimit^NumIndices;

% initialize multi-index
MultiIndex = 1;

% reverse order of the index array for easier processing
IndicesFlipped = flip( aIndices );

for iInd = 1:NumIndices
    MultiIndex = MultiIndex + ( IndicesFlipped( iInd ) - 1 ) * aLimit ^ ( iInd - 1 );
end


%% return MultiIndex, MaxMultiIndex 
end

