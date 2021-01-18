#ifndef MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP

#include "cl_GEN_Field.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete_Interpolation : virtual public Field
        {
        private:
            mtk::Mesh* mMesh;

        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Discrete_Interpolation(mtk::Mesh* aMesh);

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aNodeIndex Node index
             * @return Field value
             */
            virtual real get_field_value(uint aNodeIndex) = 0;

            /**
             * Given a node index or coordinate, returns a matrix all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index, returns a matrix of all sensitivities.
             *
             * @param aNodeIndex Node index
             * @return Matrix of sensitivities
             */
            virtual const Matrix<DDRMat>& get_field_sensitivities(uint aNodeIndex) = 0;

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations, including child nodes.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(
                    uint aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations for non-child nodes.
             *
             * @param aNodeIndex Node index
             * @return Determining ADV IDs at this node
             */
            virtual Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

        };
    }
}


#endif //MORIS_CL_GEN_FIELD_DISCRETE_INTEGRATION_HPP
