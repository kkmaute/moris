/*
 * cl_GE_Nodal_Info.hpp
 *
 *  Created on: May 17, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_PDV_INFO_HPP_
#define PROJECTS_GEN_SRC_CL_GE_PDV_INFO_HPP_

#include "catch.hpp"
#include "cl_Logger.hpp"

// GE includes
#include "cl_MTK_Mapper.hpp"
#include "cl_GE_Intersection_Object.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"

namespace moris
{
namespace ge
{
/*
 * @brief this class contains all the relevant node information associated with a {geom_rep,mesh} pair
 */
    class PDV_Info
    {
    public:
        /*
         * The constructor will require a {geom_rep,mesh} pair and will create the initial information table
         * which can then be adjusted (e.g. if additional vertices are added).
         */
        PDV_Info(std::shared_ptr< Geometry > & aGeomPointer,
                                  moris_index  aMyMeshIndex = 0 )
        {
            mMyMeshIndex = aMyMeshIndex;
            mMyGeomRep   = aGeomPointer;

            this->initialize_data_tables( aGeomPointer );
        };

        ~PDV_Info()
        {

        };
        //******************************* get functions ********************************
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
            if( mMyGeomRep->get_geom_type() == GeomType::ANALYTIC )
            {
                MORIS_ASSERT( mMyGeomRep->check_if_function_is_set(),"ge::Nodal_Info::get_field_vals() - analytic function not set " );
            }

            auto tRowItt = mMyMap.find(aNodeIndex);
            MORIS_ASSERT( tRowItt != mMyMap.end(), "ge::Nodal_Info::get_field_vals() - requested index not valid " );
            return mMyNodalFieldVals.get_row( tRowItt->second );
        };

        /*
         * @brief returns all the nodal field values
         */
        Matrix< DDRMat > get_field_vals(  )
        {
            uint tNumOfIPNodes = this->get_my_geom_rep()->get_my_mesh()->get_integration_mesh( mMyMeshIndex )->get_num_nodes();

            Matrix< DDRMat > tLSVals(tNumOfIPNodes,1, 0.0);
            for( uint n=0; n<tNumOfIPNodes; n++ )
            {
                tLSVals(n,0) = this->get_field_vals( n )( 0 );
            }
            return tLSVals;
        };

        //------------------------------------------------------------------------------
        /*
         * @brief returns the nodal sensitivities of the specified node index
         *
         *@param[in] aNodeIndex - index to the node whose information you want
         *
         *@param[out] requested information
         */
        Matrix< DDRMat > get_sensitivity_vals( moris_index aNodeIndex )
        {
            if( mMyGeomRep->get_geom_type() == GeomType::ANALYTIC )
            {
                MORIS_ASSERT( mMyGeomRep->check_if_sensitivity_function_is_set(),"ge::Nodal_Info::get_field_vals() - analytic sensitivity function not set " );
            }

            auto tCellItt = mMyMap.find(aNodeIndex);
            MORIS_ASSERT( tCellItt != mMyMap.end(), "ge::Nodal_Info::get_field_vals() - requested index not valid " );
            return mMyNodalSensitivities( tCellItt->second );
        };

        /*
         * @brief returns all the nodal sensitivity values, dphi/dp
         */
        Cell< Matrix< DDRMat > > get_sensitivity_vals(  )
        {
            uint tNumOfIPNodes = this->get_my_geom_rep()->get_my_mesh()->get_integration_mesh( mMyMeshIndex )->get_num_nodes();

            Cell< Matrix< DDRMat > > tSensitivities( tNumOfIPNodes );
            for( uint n=0; n<tNumOfIPNodes; n++ )
            {
                tSensitivities(n) = this->get_sensitivity_vals( n );
            }
            return tSensitivities;
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
            auto tCellItt = mMyMap.find(aNodeIndex);
            MORIS_ASSERT( tCellItt != mMyMap.end(), "ge::Nodal_Info::get_field_vals() - requested index not valid " );
            return mMyNodalNormals( tCellItt->second );
        };
        //------------------------------------------------------------------------------


        //------------------------------------------------------------------------------
        /*
         * @brief add node and values to the tables
         *
         */
        void add_vertex_and_value( mtk::Vertex & aVertex )
        {
            MORIS_ASSERT( mMyGeomRep->get_geom_type() == GeomType::ANALYTIC, "ge::Nodal_Info::add_vertex_and_value(): currently only set up for analytic geometry type " );

            //fixme tNewInd == end of map
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
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
                MORIS_LOG_ERROR( "note: add_vertex_and_value() - currently not set for the discrete geometry type                " );
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
                break;
            }
            case(GeomType::SDF):
            {
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
                MORIS_LOG_ERROR( "note: add_vertex_and_value() - currently not set for the SDF geometry type                     " );
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
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
        moris_index compute_intersection( Intersection_Object* aIntersectionObject )
        {
            return aIntersectionObject->compute_intersection(  );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns the local intersection point of a geometry representation with a geometric primitive
         */
        Matrix< DDRMat > get_intersection_point_local_coord( Intersection_Object* aIntersectionObject,
                                                             moris_index          aWhichIntersection = 0 )
        {
            return aIntersectionObject->get_intersection_point_local_coord( aWhichIntersection );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns the global intersection point of a geometry representation with a geometric primitive
         */
        Matrix< DDRMat > get_intersection_point_global_coord( Intersection_Object* aIntersectionObject,
                                                              moris_index          aWhichIntersection = 0 )
        {
            return aIntersectionObject->get_intersection_point_global_coord( aWhichIntersection );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief determines the intersection sensitivity of a geometry representation with a geometric primitive
         */
        void compute_intersection_sensitivity( Intersection_Object* aIntersectionObject,
                                               moris_index          aWhichIntersection = 0 )
        {
//            MORIS_ASSERT(  );
            aIntersectionObject->compute_intersection_sensitivity( aWhichIntersection );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns the intersection point sensitivity of a geometry representation with a field, dxgamma/dphi
         */
        Matrix< DDRMat > get_intersection_sensitivity( Intersection_Object* aIntersectionObject,
                                                       moris_index          aWhichIntersection = 0 )
        {
            MORIS_ASSERT( aIntersectionObject->get_field_sens_vals_cell().size() > (uint)aWhichIntersection, "ge::Nodal_Info::get_intersection_sensitivity() - sensitivity values have not been calculated or the intersection index is out of bounds " );
            return aIntersectionObject->get_field_sensitivity_vals( aWhichIntersection );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief compute and return the intersection sensitivity wrt the pdv, dxgamma/dp = dxgamma/dphi * dphi/dp
         */
        Matrix< DDRMat > get_dxgamma_dp( Intersection_Object*  aIntersectionObject,
                                         moris_index aWhichIntersection = 0 )
        {
            //fixme need an assert here to make sure the multiplication will work (dimension check)
            Matrix<DDRMat> tdxgamma_dphi = this->get_intersection_sensitivity( aIntersectionObject, aWhichIntersection );
            Matrix<DDRMat> tdphi_dp      = this->get_sensitivity_vals(  )( aWhichIntersection );
print(tdxgamma_dphi, "tdxgamma_dphi");
print(tdphi_dp, "tdphi_dp");
            return tdxgamma_dphi*tdphi_dp;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief determine LS value at a point
         */
        real compute_value_at_point( Matrix< DDRMat > aPoint )
        {
            return mMyGeomRep->get_field_val_at_coordinate( aPoint );
        }
        //------------------------------------------------------------------------------
        std::shared_ptr< Geometry > get_my_geom_rep()
        {
            return mMyGeomRep;
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
                mMyMap[n] = aGeomPointer->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mtk_vertex(n).get_index();
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
                if( mMyGeomRep->check_if_function_is_set() )
                {
                    for (uint n=0; n<aNumNodes; ++n)
                    {
                        mMyNodalFieldVals( n,0 ) = mMyGeomRep->get_field_val_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mtk_vertex(n).get_coords() );
                    }
                }
                break;
            }
            case(GeomType::DISCRETE) :
            {
                MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mesh_type() == MeshType::HMR, "Nodal_Info::create_node_val_table(): mapper for DISCRETE geom type is only set for an hmr mesh right now" );
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
                if( mMyGeomRep->check_if_sensitivity_function_is_set() )
                {
                    for (uint n=0; n<aNumNodes; ++n)
                    {
                        mMyNodalSensitivities( n ) = mMyGeomRep->get_sensitivity_dphi_dp_at_coordinate( mMyGeomRep->get_my_mesh()->get_interpolation_mesh(0)->get_mtk_vertex(n).get_coords() );
                    }
                }
                break;
            }
            case(GeomType::DISCRETE) :
            {
                /*
                 * need to be able to use finite differencing to compute derivatives for discrete class
                 */
                MORIS_LOG_ERROR( "----------------------------------------------------------------------------------------------------" );
                MORIS_LOG_ERROR( "note: sensitivity vals are currently not implemented for DISCRETE geom type, default values to zeros" );
                MORIS_LOG_ERROR( "----------------------------------------------------------------------------------------------------" );
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
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
                MORIS_LOG_ERROR( "note: sensitivity vals are currently not implemented for SDF geom type, default values to zeros" );
                MORIS_LOG_ERROR( "-----------------------------------------------------------------------------------------------" );
                Matrix< DDRMat > tZeros(mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_spatial_dim(), mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_spatial_dim(), 0.0);
                for (uint n=0; n<aNumNodes; ++n)
                {
                    mMyNodalSensitivities(n) = tZeros;
                }

                MORIS_ASSERT( mMyGeomRep->get_my_mesh()->get_interpolation_mesh( mMyMeshIndex )->get_mesh_type() == MeshType::HMR, "Nodal_Info::create_node_sensitivity_tables():  currently only set for hmr meshes" );
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
         *        * is the implementation dependent on geom_rep type?
         */
        //fixme needs to be implemented
        void create_node_normals_table()
        {
            MORIS_ASSERT(false, "ge::Nodal_Info::create_nodal_normals_table(): not implemented...yet ");
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

#endif /* PROJECTS_GEN_SRC_CL_GE_PDV_INFO_HPP_ */
