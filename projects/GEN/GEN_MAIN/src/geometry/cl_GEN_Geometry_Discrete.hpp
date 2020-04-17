//
// Created by christopherson on 4/16/20.
//

#ifndef MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
#define MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Discrete
        {
        protected:
            Cell<real*> mGeometryVariables;
            Matrix<DDRMat> mConstantParameters;

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Geometry_Discrete(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this geometry
             */
            Geometry_Discrete(Matrix<DDRMat> aConstantParameters);

        public:
            /**
             * Destructor
             */
            ~Geometry_Discrete()
            {
            }

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aEntityIndex the index of the field value
             * @return field value at the specified index
             */
            virtual real evaluate_field_value(moris_index aEntityIndex) = 0;

            /**
             * Given an index, the discrete geometry needs to return sensitivites with respect to the field value
             *
             * @param aEntityIndex the index of the field value
             * @return matrix of sensitivities
             */
            virtual Matrix<DDRMat> evaluate_sensitivity(moris_index aEntityIndex) = 0;

//        virtual GEN_Field* get_field_pointer()
//        {
//            MORIS_ASSERT(false, "GEN_Geometry::get_field_pointer() - not implemented");
//            return new GEN_Field();
//        }
///*
// * Given a node index and rank, provide the field value or minimum distance to a geometry feature
// * NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
// */
//        virtual moris::real access_field_value_with_entity_index(
//                moris::moris_index     aEntityIndex,
//                enum moris::EntityRank aEntityRank) const
//        {
//            MORIS_ASSERT(false, " access_field_value_with_entity_index not implemented. This could be due to a geometry not being based on a mesh.");
//            return 0;
//        }
//
//
///*
// * Given a local coordinate, node index and entity rank, the function returns a matrix of relevant node coordinates
// *  NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
// */
//        virtual
//        moris::Matrix< moris::DDRMat >
//        evaluate_sensitivity_dphi_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum moris::EntityRank aEntityRank)
//        {
////        std::cout<<" evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on a mesh.";
//            return moris::Matrix< moris::DDRMat >(1,1,0);
//        }
//
///*
// * Returns the advs that the provided node indices are dependent on
// */
//        virtual moris::Matrix< moris::IndexMat > get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
//        {
//            std::cout<<" get_node_adv_indices not implemented, This could be due to a geometry not being based on a mesh.";
//            return moris::Matrix< moris::IndexMat >(1,1,0);
//        }
        };
    }
}

#endif //MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
