#ifndef MORIS_CL_GEN_PROPERTY_HPP_
#define MORIS_CL_GEN_PROPERTY_HPP_

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Property : public Field
        {

        protected:
            Cell<std::shared_ptr<Property>> mPropertyDependencies;

        public:

            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aPropertyDependencies Other created properties that this property depends on
             */
            Property(Matrix<DDRMat>& aADVs,
                     Matrix<DDUMat> aPropertyVariableIndices,
                     Matrix<DDUMat> aADVIndices,
                     Matrix<DDRMat> aConstantParameters,
                     Cell<std::shared_ptr<Property>> aPropertyDependencies = Cell<std::shared_ptr<Property>>(0));

        };

    }   // end ge namespace
}   // end moris namespace

#endif /* MORIS_CL_Property_HPP_ */
