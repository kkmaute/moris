#include "cl_GEN_Single_Grain.hpp"
#include "cl_GEN_Voxel_Input.hpp"

#include "cl_Ascii.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Single_Grain::Single_Grain(
                Matrix<DDRMat>            aConstants,
                std::shared_ptr<Geometry> aVoxelGeometry,
                uint                      aIndex,
                Geometry_Field_Parameters aParameters)
                : Field(aConstants, aParameters)
                , Geometry(aParameters)
                , mVoxelGeometry(aVoxelGeometry)
                , mIndex(aIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Single_Grain::get_field_value( const Matrix<DDRMat>& aCoordinates )
        {
            moris::real tGeoValue = reinterpret_cast< Voxel_Input* >(mVoxelGeometry.get())->get_field_value(aCoordinates);

            moris::real tValue = -1.0;

            if( tGeoValue == ( real ) mIndex )
            {
                tValue =  1.0;
            }

            return tValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Single_Grain::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
        {
            MORIS_ERROR( false, "Voxel_Input::get_dfield_dadvs(), Sensitivities cannot be calculated for Voxel field.");
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
