#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(Property_Parameters aParameters)
                : mParameters(aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(std::shared_ptr<Property> aProperty)
                : mParameters(aProperty->mParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        PDV_Type Property::get_pdv_type()
        {
            return mParameters.mPDVType;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Property::is_interpolation_pdv()
        {
            return mParameters.mInterpolationPDV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDUMat> Property::get_pdv_mesh_set_indices()
        {
            return mParameters.mPDVMeshSetIndices;
        }

        //--------------------------------------------------------------------------------------------------------------

        Cell<std::string> Property::get_pdv_mesh_set_names()
        {
            return mParameters.mPDVMeshSetNames;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
