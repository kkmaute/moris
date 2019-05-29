/*
 * cl_GE_Nodal_Info.hpp
 *
 *  Created on: May 17, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_NODAL_INFO_HPP_
#define PROJECTS_GEN_SRC_CL_GE_NODAL_INFO_HPP_

#include "catch.hpp"

// GE includes
#include "cl_MTK_Mapper.hpp"
#include "cl_GE_Intersection_Object.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"

namespace moris
{
namespace ge
{
/*
 * @brief this class contains all the relevant nodal information associated with a {geom_rep,mesh} pair
 */
    class Nodal_Info
    {
    public:
        /*
         * The constructor will require a {geom_rep,mesh} pair and will create the initial information table
         * which can then be adjusted (e.g. if additional vertices are added).
         */
        Nodal_Info(std::shared_ptr< Geometry > & aGeomPointer,
                   moris_index aMyMeshIndex )
        {
            mMyMeshIndex = aMyMeshIndex;
            mMyGeomRep = aGeomPointer;

            this->initialize_data_tables( aGeomPointer );
        };

        ~Nodal_Info()
        {

        };
        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal field vals of the specified node
         *
         *@param[in] aNodeIndex - index to the node whose information you want
         *
         *@param[out] requested information
         */
        Matrix< DDRMat > get_field_vals( moris_index aNodeIndex )
        {
            auto tRowInd = mMyMap.find(aNodeIndex);
            return mMyNodalFieldVals.get_row( tRowInd->second );
        };

        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal sensitivities of the specified node
         *
         *@param[in] aNodeIndex - index to the node whose information you want
         *
         *@param[out] requested information
         */
        Matrix< DDRMat > get_sensitivity_vals( moris_index aNodeIndex )
        {
            auto tCellInd = mMyMap.find(aNodeIndex);
            return mMyNodalSensitivities( tCellInd->second );
        };

        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal normals of the specified node
         *
         *@param[in] aNodeIndex - index to the node whose information you want
         *
         *@param[out] requested information
         */
        Matrix< DDRMat > get_normal( moris_index aNodeIndex )
        {
            auto tCellInd = mMyMap.find(aNodeIndex);
            return mMyNodalNormals( tCellInd->second );
        };

        //------------------------------------------------------------------------------
        /*
         * @brief add node and values
         *
         */
        void add_vertex_and_value( mtk::Vertex & aVertex )
        {
            MORIS_ASSERT( mMyGeomRep->get_geom_type() == GeomType::ANALYTIC, "ge::Nodal_Info::add_vertex_and_value(): currently only set up for analytic geometry type" );

            moris_index tNewInd = aVertex.get_index();
            uint tSize = mMyNodalSensitivities.size();

            mMyMap[tNewInd] = tSize;

            mMyNodalFieldVals.resize(tSize+1,1);  // need to know number of fields to set this size, assuming only one for now

            switch(mMyGeomRep->get_geom_type())
            {
            case(GeomType::ANALYTIC):
            {
                mMyNodalFieldVals(tSize,0) = mMyGeomRep->get_field_val_at_coordinate( aVertex.get_coords() );
    //            mMyNodalFieldVals(tSize,0) = mMyGeomRep->get_field_val_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(tNewInd).get_coords() );
                mMyNodalSensitivities.push_back(mMyGeomRep->get_sensitivity_dphi_dp_at_coordinate( aVertex.get_coords() ));
    //            mMyNodalSensitivities.push_back(mMyGeomRep->get_sensitivity_dphi_dp_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(tNewInd).get_coords() ));
                // have not implemented node normals yet.....
    //            mMyNodalNormals.pushback( ) = ...
                break;
            }
            case(GeomType::DISCRETE):
            {
                break;
            }
            case(GeomType::SDF):
            {
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Nodal_Info::add_vertex_and_value(): geometry type not supported " );
                break;
            }

            }
        };
        //------------------------------------------------------------------------------
        /*
         * @brief determines the intersection point of a geometry representation with a geometric primitive
         */
        void compute_intersection( Intersection_Object_Line* aIntersectionObject )
        {
            MORIS_ASSERT( aIntersectionObject->get_my_global_coord().n_cols() != 0, "ge::Nodal_Info::determine_intersection(): global coordinates of Intersection_Object not set" );
            MORIS_ASSERT( aIntersectionObject->get_my_time_coord().n_cols() != 0, "ge::Nodal_Info::determine_intersection(): time coordinates of Intersection_Object not set" );
            MORIS_ASSERT( aIntersectionObject->get_my_field_vals().n_cols() != 0, "ge::Nodal_Info::determine_intersection(): field values of Intersection_Object not set" );

            Matrix< DDRMat > tParamPoint = aIntersectionObject->get_my_param_point();
            for(uint k=0; k<201; ++k)
            {
                if (k == 200)
                {
                    std::cout<<"no intersection point found"<<std::endl;
                    break;
                }
                aIntersectionObject->get_my_field_interp()->set_space_time( tParamPoint );

                Matrix< DDRMat > tRes = aIntersectionObject->get_my_field_interp()->val();     // R = phi - 0
                if (tRes(0,0) == Approx(0.0))
                {
                    aIntersectionObject->is_intersected();
                    aIntersectionObject->set_intersection_point( aIntersectionObject->get_my_geom_interp()->valx(tParamPoint) );
std::cout<<"-2-2-2-2-2-2-2-2-2-2-2-2-22-2-2-2-2-2-2-2"<<std::endl;
                    break;
                }
                tParamPoint(0) += 0.01;
            };
        }
        //------------------------------------------------------------------------------
        /*
         * @brief determine LS value at a point
         */
        real compute_value_at_point( Matrix< DDRMat > aPoint )
        {
            return mMyGeomRep->get_field_val_at_coordinate( aPoint );
        }
    private:
        //------------------------------------------------------------------------------
        /*
         * @brief create data tables and relationship between nodes and data
         *
         */
        void initialize_data_tables( std::shared_ptr< Geometry > & aGeomPointer )
        {
            uint tNumNodes = aGeomPointer->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_num_nodes();
            // build relations map
            for(uint n=0; n<tNumNodes; ++n)
            {
                mMyMap[n] = aGeomPointer->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(n).get_index();
            }

            //fixme need to check how many fields are on the mesh and initialize that many columns for the mMyNodalFieldVals matrix
            mMyNodalFieldVals.set_size(tNumNodes,1);
            mMyNodalSensitivities.resize(tNumNodes);

            this->create_node_val_table( tNumNodes );
            this->create_node_sensitivity_table( tNumNodes );
            // compute nodal normals currently not implemented
//            this->create_nodal_normals_table();
        };

        //------------------------------------------------------------------------------
        /*
         * @brief compute the field vals at the nodes and create the field val table
         *        * compute field vals at nodes - store in table, register table
         *        * implementation is dependent on geom_rep type
         */
        void create_node_val_table( uint aNumNodes )
        {
            switch(mMyGeomRep->get_geom_type())
            {

            case(GeomType::ANALYTIC) :
            {
                mMyGeomRep->check_if_functions_are_set();
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalFieldVals( n,0 ) = mMyGeomRep->get_field_val_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(n).get_coords() );
                }
                break;
            }
            case(GeomType::DISCRETE) :
            {
                MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(mMyMeshIndex)->get_mesh_type() == MeshType::HMR, "Nodal_Info::create_node_val_table(): mapper for DISCRETE geom type is only set for an hmr mesh right now" );
                /*
                 * terminology is kind of confusing:
                 *  *mMyGeomRep->get_target_field() gives the target field to map
                 *  *mMyGeomRep->get_output_field() gives the output field resulting from the map
                 */
                mapper::Mapper tMyMapper( mMyGeomRep->get_my_mesh(), mMyMeshIndex, mMyGeomRep->get_my_target_field()->get_bspline_order() );
                // default is to map from a node field to a B-spline field, this can be changed if necessary
                tMyMapper.perform_mapping( mMyGeomRep->get_my_target_field()->get_label(), EntityRank::NODE, mMyGeomRep->get_my_output_field()->get_label(), mMyGeomRep->get_my_output_field()->get_bspline_rank() );

                mMyNodalFieldVals = mMyGeomRep->get_my_output_field()->get_coefficients();
                break;
            }
            case(GeomType::SDF) :
            {
                MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "Nodal_Info::create_node_val_tables():  currently only set for hmr meshes" );

                for(uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalFieldVals(n,0) = mMyGeomRep->get_sdf_vals(n);
                }
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Nodal_Info::create_node_val_table(): geometry type not supported " );
            }

            }
        };

        //------------------------------------------------------------------------------
        /*
         * @brief compute the field sensitivities at the nodes and create the field sensitivity table
         *        * compute field sensitivities at nodes - store in table, register table
         *        * implementation is dependent on geom_rep type
         */
        void create_node_sensitivity_table( uint aNumNodes )
        {
            switch(mMyGeomRep->get_geom_type())
            {
            case(GeomType::ANALYTIC) :
            {
                mMyGeomRep->check_if_functions_are_set();
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalSensitivities( n ) = mMyGeomRep->get_sensitivity_dphi_dp_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(n).get_coords() );
                }
                break;
            }
            case(GeomType::DISCRETE) :
            {
                /*
                 * need to be able to use finite differencing to compute derivatives for discrete class
                 */
                std::cout<<"----------------------------------------------------------------------------------------------------"<<std::endl;
                std::cout<<"note: sensitivity vals are currently not implemented for DISCRETE geom type, default values to zeros"<<std::endl;
                std::cout<<"----------------------------------------------------------------------------------------------------"<<std::endl;
                Matrix< DDRMat > tZeros(mMyGeomRep->get_my_mesh()->get_interpolation_mesh(mMyMeshIndex)->get_spatial_dim(), mMyGeomRep->get_my_mesh()->get_interpolation_mesh(mMyMeshIndex)->get_spatial_dim(), 0.0);
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalSensitivities(n) = tZeros;
                }
                break;
            }
            case(GeomType::SDF) :
            {
                /*
                 * not sure how getting the sensitivities from the SDF type will work (via raycast?)
                 */
                std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
                std::cout<<"note: sensitivity vals are currently not implemented for SDF geom type, default values to zeros"<<std::endl;
                std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
                Matrix< DDRMat > tZeros(mMyGeomRep->get_my_mesh()->get_interpolation_mesh(mMyMeshIndex)->get_spatial_dim(), mMyGeomRep->get_my_mesh()->get_interpolation_mesh(mMyMeshIndex)->get_spatial_dim(), 0.0);
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalSensitivities(n) = tZeros;
                }

                MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR, "Nodal_Info::create_node_sensitivity_tables():  currently only set for hmr meshes" );
                break;
            }
            default :
            {
                MORIS_ASSERT( false, "Nodal_Info::create_node_sensitivity_table(): geometry type not supported " );
            }

            }
        };

        //------------------------------------------------------------------------------
        /*
         * @brief compute the nodal normals and create table
         *        * compute normals at nodes (geometric and facet) - store in table, register table
         *        * is the implementation is dependent on geom_rep type?
         */
        //fixme needs to be implemented
        void create_node_normals_table()
        {
            MORIS_ASSERT(false, "ge::Nodal_Info::create_nodal_normals_table(): not yet implemented");
        };
        //------------------------------------------------------------------------------
        std::shared_ptr< Geometry > mMyGeomRep;
        moris_index                 mMyMeshIndex;

        std::unordered_map< moris_index, moris_index > mMyMap; // (vertex index, row/cell index for data tables)

        Matrix< DDRMat >                mMyNodalFieldVals;
        moris::Cell< Matrix< DDRMat > > mMyNodalSensitivities;
        moris::Cell< Matrix< DDRMat > > mMyNodalNormals;
        //------------------------------------------------------------------------------
    protected:

    };
} /* namespace gen */
} /* namespace moris */

#endif /* PROJECTS_GEN_SRC_CL_GE_NODAL_INFO_HPP_ */
