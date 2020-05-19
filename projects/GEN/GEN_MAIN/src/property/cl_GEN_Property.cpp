#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        Property::Property(Matrix<DDRMat>& aADVs,
                           Matrix<DDUMat> aPropertyVariableIndices,
                           Matrix<DDUMat> aADVIndices,
                           Matrix<DDRMat> aConstantParameters,
                           Cell<std::shared_ptr<Property>> aPropertyDependencies)
                           : Field(aADVs, aPropertyVariableIndices, aADVIndices, aConstantParameters),
                           mPropertyDependencies(aPropertyDependencies)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }   // end ge namespace
}   // end moris namespace
