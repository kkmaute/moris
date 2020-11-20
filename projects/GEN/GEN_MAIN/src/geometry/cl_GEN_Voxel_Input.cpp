#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Voxel_Input::Voxel_Input(
                Matrix<DDRMat>&                aADVs,
                Matrix<DDUMat>                 aGeometryVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstantParameters,
                std::string                    aVoxelFieldName,
                std::string                    aName,
                Matrix<DDSMat>                 aNumRefinements,
                Matrix<DDSMat>                 aRefinementMeshIndices,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex)
        : Field(aADVs,
                aGeometryVariableIndices,
                aADVIndices,
                aConstantParameters,
                aName,
                aNumRefinements,
                aRefinementMeshIndices,
                aRefinementFunctionIndex,
                aBSplineMeshIndex,
                0.0,
                0.0)
        {
            this->read_voxel_data( aVoxelFieldName );
        }

        //--------------------------------------------------------------------------------------------------------------

        Voxel_Input::Voxel_Input(
                sol::Dist_Vector*              aOwnedADVs,
                Matrix<DDUMat>                 aGeometryVariableIndices,
                Matrix<DDUMat>                 aADVIndices,
                Matrix<DDRMat>                 aConstantParameters,
                std::string                    aVoxelFieldName,
                std::string                    aName,
                Matrix<DDSMat>                 aNumRefinements,
                Matrix<DDSMat>                 aRefinementMeshIndices,
                sint                           aRefinementFunctionIndex,
                sint                           aBSplineMeshIndex)
        : Field(aOwnedADVs,
                aGeometryVariableIndices,
                aADVIndices,
                aConstantParameters,
                aName,
                aNumRefinements,
                aRefinementMeshIndices,
                aRefinementFunctionIndex,
                aBSplineMeshIndex,
                0.0,
                0.0)
        {
            this->read_voxel_data( aVoxelFieldName );
        }

        //--------------------------------------------------------------------------------------------------------------

        Voxel_Input::Voxel_Input(
                Matrix<DDRMat>           aConstantParameters,
                std::string                    aVoxelFieldName,
                std::string              aName,
                Matrix<DDSMat>           aNumRefinements,
                Matrix<DDSMat>           aRefinementMeshIndices,
                sint                     aRefinementFunctionIndex,
                sint                     aBSplineMeshIndex)
        : Field(aConstantParameters,
                aName,
                aNumRefinements,
                aRefinementMeshIndices,
                aRefinementFunctionIndex,
                aBSplineMeshIndex,
                0.0,
                0.0)
        {
            this->read_voxel_data( aVoxelFieldName );
        }

        //--------------------------------------------------------------------------------------------------------------

        real Voxel_Input::get_field_value( const Matrix<DDRMat>& aCoordinates )
        {
            moris::Matrix< DDRMat > tDomainDimensions = { {4.0}, {4.0}, {4.0} };
            moris::Matrix< DDRMat > tDomainOffset = { {-2.0}, {-2.0}, {-2.0} };

            moris::real tVoxelSizeX = tDomainDimensions( 0 ) / mVoxelsInX;
            moris::real tVoxelSizeY = tDomainDimensions( 1 ) / mVoxelsInY;
            moris::real tVoxelSizeZ = tDomainDimensions( 2 ) / mVoxelsInZ;

            moris::uint tI = std::floor( ( aCoordinates( 0 ) - tDomainOffset( 0 ) ) / tVoxelSizeX );
            moris::uint tJ = std::floor( ( aCoordinates( 1 ) - tDomainOffset( 1 ) ) / tVoxelSizeY );
            moris::uint tK = std::floor( ( aCoordinates( 2 ) - tDomainOffset( 2 ) ) / tVoxelSizeZ );

            if( aCoordinates( 0 ) >=  tDomainDimensions( 0 ) + tDomainOffset( 0  ))
            {
                tI = mVoxelsInX-1;
            }
            if( aCoordinates( 1 ) >=  tDomainDimensions( 1 ) + tDomainOffset( 1 ))
            {
                tJ = mVoxelsInY-1;
            }
            if( aCoordinates( 2 ) >=  tDomainDimensions( 2 ) + tDomainOffset( 2 ))
            {
                tK = mVoxelsInZ-1;
            }

            // mVoxelField columns are ordered - VoxelIndex - GainsId - I - J - K
            moris::uint tRow =
                    tI * mVoxelsInY * mVoxelsInZ +
                    tJ * mVoxelsInZ +
                    tK;

            return ( moris::real )mVoxelField( tRow, 1 );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Voxel_Input::get_field_sensitivities(const Matrix<DDRMat>& aCoordinates)
        {
            MORIS_ERROR( false, "Voxel_Input::get_field_sensitivities(), Sensitivities cannot be calculated for Voxel field.");
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Voxel_Input::read_voxel_data( std::string aVoxelFieldName)
        {
            // build Ascii reader
            moris::Ascii tAsciiReader( aVoxelFieldName, moris::FileMode::OPEN_RDONLY );

            // get number of lines in asci file
            moris::uint tNumLines = tAsciiReader.length();

            // set matrix for voxel field
            mVoxelField.set_size( tNumLines, 5, MORIS_UINT_MAX );

            for( uint Ik = 0; Ik < tNumLines; Ik++ )
            {
                std::string tFileLine = tAsciiReader.line( Ik );

                // reset position
                size_t tPos = 0;

                if( !tFileLine.empty()  )
                {
                    uint tCountCol = std::count( tFileLine.begin(), tFileLine.end(), ' ') + 1;

                    for( uint Ii = 0; Ii < tCountCol-1; Ii++ )
                    {
                        // find string
                        tPos = tFileLine.find( " " );

                        // copy value into output matrix
                        if( tPos < tFileLine.size() )
                        {
                            mVoxelField( Ik, Ii ) = stod(  tFileLine.substr( 0, tPos ) );
                            tFileLine =  tFileLine.substr( tPos+1, tFileLine.size() );
                        }
                    }

                    // copy value into output matrix
                    if( !tFileLine.empty() )
                    {
                        mVoxelField( Ik, tCountCol-1 ) = stod(  tFileLine.substr( 0, tFileLine.size() ) );
                    }
                }
            }

            MORIS_ERROR( mVoxelField.max() != MORIS_UINT_MAX, "Voxel_Input::Voxel_Input(), Matrix not correctly initialized");

            mNumGrainInd = mVoxelField.get_column( 1 ).max();

            // Voxel indices start with 1
            mVoxelsInX   = mVoxelField.get_column( 2 ).max() + 1;
            mVoxelsInY   = mVoxelField.get_column( 3 ).max() + 1;
            mVoxelsInZ   = mVoxelField.get_column( 4 ).max() + 1;
        }


        //--------------------------------------------------------------------------------------------------------------

    }
}
