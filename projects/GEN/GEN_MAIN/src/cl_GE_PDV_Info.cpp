/*
 * cl_GE_Factory.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */
#include "cl_GE_PDV_Info.hpp"

#include "cl_MTK_Mapper.hpp"
#include "cl_GE_Geometry.hpp"


namespace moris
{
    namespace ge
    {
    void PDV_Info::create_node_val_table( uint aNumNodes, moris_index aSubIndex )
    {
        switch(mMyGeomRep->get_geom_type())
        {

        case(GeomType::ANALYTIC) :
        {
            if( mMyGeomRep->check_if_function_is_set( aSubIndex ) )
            {
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalFieldVals( n,aSubIndex ) = mMyGeomRep->get_field_val_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mtk_vertex(n).get_coords(), aSubIndex );
                }
            }
            break;
        }
        case(GeomType::DISCRETE) :
        {
            MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mesh_type() == MeshType::HMR, "PDV_Info::create_node_val_table() - mapper for DISCRETE geom type is only set for an hmr mesh right now" );
            /*
             * terminology is kind of confusing:
             *  *mMyGeomRep->get_target_field() gives the target field to map
             *  *mMyGeomRep->get_output_field() gives the output field resulting from the map
             */
            mapper::Mapper tMyMapper( mMyGeomRep->get_my_mesh(), mMyMeshIndex, mMyGeomRep->get_my_target_field()->get_bspline_order() );
            // default is to map from a node field to a B-spline field, this can be changed if necessary
            tMyMapper.perform_mapping( mMyGeomRep->get_my_target_field()->get_label(), EntityRank::NODE, mMyGeomRep->get_my_output_field()->get_label(), mMyGeomRep->get_my_output_field()->get_bspline_rank() );

            //fixme: get_coefficients() returns a matrix of values, need to incorporate these into the node value table somehow...
            mMyNodalFieldVals = mMyGeomRep->get_my_output_field()->get_coefficients();
            break;
        }
        case(GeomType::SDF) :
        {
            MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "PDV_Info::create_node_val_tables() - currently only set for hmr meshes" );

            for(uint n=0; n<aNumNodes; ++n)
            {
                mMyNodalFieldVals(n,0) = mMyGeomRep->get_sdf_vals(n);
            }
            break;
        }
        default :
        {
            MORIS_ERROR( false, "PDV_Info::create_node_val_table() - geometry type not supported " );
        }

        }
    }
    } /* namespace gen */
} /* namespace moris */
