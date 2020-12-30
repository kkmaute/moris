#ifndef MORIS_CL_GEN_SCALED_FIELD_HPP
#define MORIS_CL_GEN_SCALED_FIELD_HPP

#include "cl_GEN_Property.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Scaled_Field : public Property
        {

        public:
            /**
             * Constructor
             *
             * @param aADVs Reference to the full advs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aParameters Additional parameters
             */
            Scaled_Field(
                    Matrix<DDRMat>&              aADVs,
                    Matrix<DDUMat>               aPropertyVariableIndices,
                    Matrix<DDUMat>               aADVIndices,
                    Matrix<DDRMat>               aConstants,
                    Cell<std::shared_ptr<Field>> aFieldDependencies,
                    Field_Parameters             aParameters = {});

            /**
             * Constructor
             *
             * @param aOwnedADVs Distributed owned ADVs
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aFieldDependencies Other created fields that this property depends on
             * @param aParameters Additional parameters
             */
            Scaled_Field(
                    sol::Dist_Vector*            aOwnedADVs,
                    Matrix<DDUMat>               aPropertyVariableIndices,
                    Matrix<DDUMat>               aADVIndices,
                    Matrix<DDRMat>               aConstants,
                    Cell<std::shared_ptr<Field>> aFieldDependencies,
                    Field_Parameters             aParameters = {});

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Property value
             */
            real get_field_value(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, evaluates the sensitivity of the property field with respect to all of the
             * property variables.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);
        };
    }
}

#endif //MORIS_CL_GEN_SCALED_FIELD_HPP
