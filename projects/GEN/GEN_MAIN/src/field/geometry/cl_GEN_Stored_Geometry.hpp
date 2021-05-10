#ifndef MORIS_CL_GEN_STORED_GEOMETRY_HPP
#define MORIS_CL_GEN_STORED_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace ge
    {
        class Stored_Geometry : public Geometry, public Field_Discrete_Integration
        {

        private:
            std::shared_ptr<Geometry> mGeometry;
            mtk::Mesh* mMesh;
            Matrix<DDRMat> mFieldValues;

        public:

            /**
             * Constructor
             *
             * @param aMesh The mesh pointer where node information can be obtained
             * @param aGeometry Geometry for obtaining field values to store
             */
            Stored_Geometry(
                    mtk::Mesh*                aMesh,
                    std::shared_ptr<Geometry> aGeometry);

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Field value
             */
            real get_field_value(uint aNodeIndex);

            /**
             * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_dfield_dadvs(uint aNodeIndex);

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @return Determining ADV IDs at this node
             */
            Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

            /**
             * Given a node index, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aNodeIndex Node index
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            void get_dfield_dcoordinates(
                    uint            aNodeIndex,
                    Matrix<DDRMat>& aSensitivities);

            /**
             * Resets all nodal information about geometry field values.
             */
            void reset_nodal_data();

        private:

            /**
             * Evaluates and stores the nodal values of this geometry for use later.
             */
            void evaluate_nodal_values();

        };
    }
}

#endif //MORIS_CL_GEN_STORED_GEOMETRY_HPP
