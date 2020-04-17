//
// Created by christopherson on 4/16/20.
//

#ifndef MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
#define MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP

namespace moris
{
    namespace ge
    {
        virtual GEN_Field* get_field_pointer()
        {
            MORIS_ASSERT(false, "GEN_Geometry::get_field_pointer() - not implemented");
            return new GEN_Field();
        }
/*
 * Given a node index and rank, provide the field value or minimum distance to a geometry feature
 * NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
 */
        virtual moris::real access_field_value_with_entity_index(
                moris::moris_index     aEntityIndex,
                enum moris::EntityRank aEntityRank) const
        {
            MORIS_ASSERT(false, " access_field_value_with_entity_index not implemented. This could be due to a geometry not being based on a mesh.");
            return 0;
        }


/*
 * Given a local coordinate, node index and entity rank, the function returns a matrix of relevant node coordinates
 *  NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
 */
        virtual
        moris::Matrix< moris::DDRMat >
        evaluate_sensitivity_dphi_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum moris::EntityRank aEntityRank)
        {
//        std::cout<<" evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on a mesh.";
            return moris::Matrix< moris::DDRMat >(1,1,0);
        }

/*
 * Returns the advs that the provided node indices are dependent on
 */
        virtual moris::Matrix< moris::IndexMat > get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
        {
            std::cout<<" get_node_adv_indices not implemented, This could be due to a geometry not being based on a mesh.";
            return moris::Matrix< moris::IndexMat >(1,1,0);
        }
    }
}

#endif //MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
