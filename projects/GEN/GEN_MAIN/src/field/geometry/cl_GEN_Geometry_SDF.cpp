#include "cl_GEN_Geometry_SDF.hpp"

// SDF
#include "cl_SDF_Generator.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry_SDF::Geometry_SDF( std::string tObjectPath,
                                    Geometry_Field_Parameters aParameters)
                : Field(Matrix<DDRMat>(1, 1, 0.0), aParameters)
                , Geometry(aParameters)
                , Field_Discrete_Integration()
        {
            mObjectPath = tObjectPath;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Geometry_SDF::get_field_value(uint aNodeIndex)
        {
            MORIS_ASSERT( mValues.numel() > aNodeIndex, "node index out of bounce, mValues numel: %i , aNodeIndex %i",
                    mValues.numel(),
                    aNodeIndex );
            return mValues( aNodeIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

//        const Matrix<DDRMat>& Geometry_SDF::get_dfield_dadvs(uint aNodeIndex)
//        {
//            MORIS_ERROR(false, "get_dfield_dadvs function is not implemented for a mesh field geometry.");
//            return mSensitivities;
//        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_SDF::evaluate_nodal_values()
        {
            mtk::Mesh* tMesh = mMeshPair.get_interpolation_mesh();

            sdf::SDF_Generator tSDFGenerator( mObjectPath, true );

            tSDFGenerator.calculate_sdf( tMesh, mValues );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry_SDF::reset_nodal_data()
        {
            // Reset child nodes
            Field_Discrete_Integration::reset_nodal_data();

            // Re-evaluate field values
            this->evaluate_nodal_values();
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
