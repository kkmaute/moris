#ifndef MORIS_CL_GEN_SCALED_FIELD_HPP
#define MORIS_CL_GEN_SCALED_FIELD_HPP

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Discrete_Interpolation.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Scaled_Field : public Property, public Field_Discrete_Interpolation
        {

        private:
            std::shared_ptr<Field> mField;

        public:
            /**
             * Constructor
             *
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the property variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aField Field that this property will scale
             * @param aParameters Additional parameters
             */
            template <typename Vector_Type>
            Scaled_Field(
                    Vector_Type&              aADVs,
                    Matrix<DDUMat>            aPropertyVariableIndices,
                    Matrix<DDUMat>            aADVIndices,
                    Matrix<DDRMat>            aConstants,
                    std::shared_ptr<Field>    aField,
                    Property_Field_Parameters aParameters = {})
                    : Field(aADVs, aPropertyVariableIndices, aADVIndices, aConstants, aParameters)
                    , Property(aParameters)
                    , mField(aField)
            {
                MORIS_ERROR(mFieldVariables.size() == 1, "A scaled field property must have one scaling factor.");
                MORIS_ERROR(aPropertyVariableIndices.length() == 0,
                            "A scaled field property must have a constant scaling factor for now.");
            }

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Property value
             */
            real get_base_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, evaluates the sensitivity of the property field with respect to all of the
             * property variables.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_base_dfield_dadvs(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_base_determining_adv_ids(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

        private:

            /**
             * Sets the dependencies of this property after they have been found by update_dependencies(). By default
             * does nothing.
             *
             * @param aDependencyFields Fields that this property depends on.
             */
            void set_dependencies(Cell<std::shared_ptr<Field>> aDependencyFields);

        };
    }
}

#endif //MORIS_CL_GEN_SCALED_FIELD_HPP
