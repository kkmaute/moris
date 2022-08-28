%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function save_matrix_to_binary( aMatrix, aPath )
    tFileInfo = dir(aPath);
    
    tNumberOfRows = size( aMatrix, 1 );
    tNumberOfCols = size( aMatrix, 2 );
    
    % flatten matrix to vector
    tVector = zeros( tNumberOfRows * tNumberOfCols, 1 );
    
    tCount = 0;
    for ( j = 1:tNumberOfCols )
        for( i = 1:tNumberOfRows )
           tCount = tCount + 1;
           tVector( tCount ) = aMatrix( i, j );
        end
    end
    
    tFileID=fopen( aPath, "w" );
    fwrite( tFileID, tNumberOfRows, "uint32" );
    fwrite( tFileID, tNumberOfCols, "uint32" );
    fwrite( tFileID, tVector, "double");
    fclose(tFileID);
end


