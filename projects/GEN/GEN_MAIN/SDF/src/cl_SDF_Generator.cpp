
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Core.hpp"

#include "cl_SDF_Generator.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        SDF_Generator::SDF_Generator(
                const std::string & aObjectPath,
                const bool aVerboseFlag ) :
            mObject( aObjectPath ),
            mVerboseFlag ( aVerboseFlag )
        {

        }

//-------------------------------------------------------------------------------

        SDF_Generator::SDF_Generator(
                const std::string & aObjectPath,
                Matrix<DDRMat>   aObjectOffset,
                const bool aVerboseFlag  ) :
                    mObject( aObjectPath, aObjectOffset ),
                    mVerboseFlag ( aVerboseFlag )
        {
        }

//-------------------------------------------------------------------------------

        void
        SDF_Generator::raycast(
                mtk::Mesh          * aMesh,
                Matrix< IndexMat > & aElementsAtSurface )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData );

            // perform raycast
            tCore.calculate_raycast( aElementsAtSurface );
        }
//-------------------------------------------------------------------------------

        void
        SDF_Generator::raycast(
                std::shared_ptr< mtk::Mesh > aMesh,
                Matrix< IndexMat > & aElementsAtSurface )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData );

            // perform raycast
            tCore.calculate_raycast( aElementsAtSurface );
        }

//-------------------------------------------------------------------------------

        void
        SDF_Generator::raycast(
                mtk::Mesh          * aMesh,
                Matrix< IndexMat > & aElementsAtSurface,
                Matrix< IndexMat > & aElementsInVolume )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData );

            // perform raycast
            tCore.calculate_raycast( aElementsAtSurface, aElementsInVolume );
        }

//-------------------------------------------------------------------------------

        void
        SDF_Generator::raycast(
                std::shared_ptr< mtk::Mesh > aMesh,
                Matrix< IndexMat > & aElementsAtSurface,
                Matrix< IndexMat > & aElementsInVolume )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh, mVerboseFlag );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData, mVerboseFlag );

            // perform raycast
            tCore.calculate_raycast( aElementsAtSurface, aElementsInVolume );
        }

//-------------------------------------------------------------------------------

        void
        SDF_Generator::calculate_sdf(
                mtk::Mesh          * aMesh,
                Matrix< DDRMat>    & aSDF )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh, mVerboseFlag );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData, mVerboseFlag );

            // calculate SDF
            tCore.calculate_raycast_and_sdf( aSDF );

            //tCore.save_to_vtk( "sdf_mesh.vtk");
        }
//-------------------------------------------------------------------------------

        void
        SDF_Generator::calculate_sdf(
                std::shared_ptr< mtk::Mesh > aMesh,
                Matrix< DDRMat>    & aSDF )
        {
            // create mesh wrapper
            Mesh tMesh( aMesh, mVerboseFlag );

            // create data container
            Data tData( mObject );

            // create core
            Core tCore( tMesh, tData, mVerboseFlag );

            // calculate SDF
            tCore.calculate_raycast_and_sdf( aSDF );
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
