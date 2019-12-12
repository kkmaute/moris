/*
 * cl_MTK_Element_Set.hpp
 *
 *  Created on: Jun 24, 2019
 *      Author: schmidt
 */

#ifndef SRC_MESH_CL_MTK_SET_HPP_
#define SRC_MESH_CL_MTK_SET_HPP_

#include <string>

#include "typedefs.hpp" //MRS/COR/src
#include "fn_unique.hpp" //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

#include "cl_MTK_Cell_Cluster.hpp" //MTK/src
#include "cl_MTK_Side_Cluster.hpp" //MTK/src
#include "cl_Communication_Tools.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
        class Set
        {
        private :
            // name of the set
            std::string mSetName;

            // interpolation mesh geometry type
            mtk::Geometry_Type mIPGeometryType = mtk::Geometry_Type::UNDEFINED;

            // space interpolation order for IP cells
            mtk::Interpolation_Order mIPSpaceInterpolationOrder;

            // space interpolation order for IG cells
            mtk::Interpolation_Order mIGSpaceInterpolationOrder;

            bool mIsTrivialMaster = false;
            bool mIsTrivialSlave = false;

            bool mMasterLock = false;
            bool mSlaveLock = false;

//------------------------------------------------------------------------------

        protected :

            moris::Cell< Cluster const * > mSetClusters;

            // integration mesh geometry type
            mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

//------------------------------------------------------------------------------

        private:

            moris::uint mSpatialDim;

//------------------------------------------------------------------------------

        protected:

            void communicate_ip_geometry_type()
            {
                mtk::Geometry_Type tIPGeometryType = mtk::Geometry_Type::UNDEFINED;

                if( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    tIPGeometryType = mSetClusters( 0 )->get_interpolation_cell().get_geometry_type();
                }

                uint tRecIPGeometryType = (uint) mtk::Geometry_Type::UNDEFINED;

                min_all( (uint)tIPGeometryType, tRecIPGeometryType );

                mIPGeometryType = static_cast<enum mtk::Geometry_Type> (tRecIPGeometryType);

//                MORIS_ASSERT( mIPGeometryType != mtk::Geometry_Type::UNDEFINED, " communicate_type(); undefined geometry type on all processors");
            };

//------------------------------------------------------------------------------

            // FIXME should be userdefined in FEM
            void communicate_interpolation_order()
            {
//                MORIS_ASSERT( mIPGeometryType != mtk::Geometry_Type::UNDEFINED,
//                        " communicate_interpolation_order(); undefined geometry type on this processor. Try calling communicate_ip_geometry_type() first.");

                mtk::Interpolation_Order tIPInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;
                mtk::Interpolation_Order tIGInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

                if( mSetClusters.size() > 0 )
                {
                    // interpolation order for IP cells fixme
                    tIPInterpolationOrder = mSetClusters( 0 )->get_interpolation_cell( mtk::Master_Slave::MASTER ).get_interpolation_order();

                    // interpolation order for IG cells fixme
                    tIGInterpolationOrder = mSetClusters( 0 )->get_primary_cells_in_cluster( mtk::Master_Slave::MASTER )( 0 )->get_interpolation_order();
                }

                uint tRecIPInterpolationOrder = (uint) mtk::Interpolation_Order::UNDEFINED;
                uint tRecIGInterpolationOrder = (uint) mtk::Interpolation_Order::UNDEFINED;

                min_all( (uint)tIPInterpolationOrder, tRecIPInterpolationOrder );
                min_all( (uint)tIGInterpolationOrder, tRecIGInterpolationOrder );

                mIPSpaceInterpolationOrder = static_cast<enum mtk::Interpolation_Order> (tRecIPInterpolationOrder);
                mIGSpaceInterpolationOrder = static_cast<enum mtk::Interpolation_Order> (tRecIGInterpolationOrder);

//                MORIS_ASSERT( mIPSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " communicate_interpolation_order(); undefined ip interpolation order on this processor");
//                MORIS_ASSERT( mIGSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " communicate_interpolation_order(); undefined ig interpolation order on this processor");
            }

//------------------------------------------------------------------------------

            // FIXME should be userdefined in FEM
            void communicate_is_trivial_flag( const mtk::Master_Slave aIsMaster )
            {
                sint tIsTrivial = 1;

                if( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    tIsTrivial = (sint)mSetClusters( 0 )->is_trivial( aIsMaster ); //FIXME change for double sided set
                }

                sint tIsTrivialMax = (sint) false;

                max_all( tIsTrivial, tIsTrivialMax );

                if( tIsTrivialMax == 1 && aIsMaster == mtk::Master_Slave::MASTER )
                {
                    mIsTrivialMaster = true;
                }
                else if( tIsTrivialMax == 1 && aIsMaster == mtk::Master_Slave::SLAVE )
                {
                    mIsTrivialSlave = true;
                }
            };

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Set()
            { };

            Set(std::string                   aName,
                moris::Cell<Cluster const *>  aBlockSetClusters,
                const uint                    aSpatialDim ) : mSetName( aName ),
                                                              mSetClusters( aBlockSetClusters ),
                                                              mSpatialDim( aSpatialDim )
            {
                this->communicate_ip_geometry_type();

                this->communicate_interpolation_order();
            };

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Set(){};

//------------------------------------------------------------------------------

            std::string get_set_name()
            {
                return mSetName;
            }

//------------------------------------------------------------------------------

            uint get_clusters() const
            {
                MORIS_ERROR( !(mSpatialDim < 1) || !(mSpatialDim > 3), "Set::get_spatial_dim(), Spatial dim < 1 or > 3" );
                return mSpatialDim;
            }

//------------------------------------------------------------------------------

            uint get_spatial_dim() const
            {
                MORIS_ERROR( !(mSpatialDim < 1) || !(mSpatialDim > 3), "Set::get_spatial_dim(), Spatial dim < 1 or > 3" );
                return mSpatialDim;
            }

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
//              virtual const moris::Matrix< DDUMat > &
//              get_list_of_block_cell_clusters() const = 0;

            bool is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                if ( !mMasterLock && aIsMaster == mtk::Master_Slave::MASTER )
                {
                    this->communicate_is_trivial_flag( aIsMaster );

                    mMasterLock = true;
                }
                else if ( !mSlaveLock && aIsMaster == mtk::Master_Slave::SLAVE )
                {
                    this->communicate_is_trivial_flag( aIsMaster );

                    mSlaveLock = true;
                }

                if ( aIsMaster == mtk::Master_Slave::MASTER )
                {
                    return mIsTrivialMaster;
                }
                else if ( aIsMaster == mtk::Master_Slave::SLAVE )
                {
                    return mIsTrivialSlave;
                }
                else
                {
                    MORIS_ASSERT( false, " is_trivial(); undefined type. has to be master or slave");

                    return false;
                }
            };

//------------------------------------------------------------------------------

            mtk::Geometry_Type get_interpolation_cell_geometry_type()
            {
                return mIPGeometryType;
            }

//------------------------------------------------------------------------------

            mtk::Geometry_Type get_integration_cell_geometry_type()
            {
                return mIGGeometryType;
            }

//------------------------------------------------------------------------------

            mtk::Interpolation_Order get_interpolation_cell_interpolation_order()
            {
                return mIPSpaceInterpolationOrder;
            }

//------------------------------------------------------------------------------

            mtk::Interpolation_Order get_integration_cell_interpolation_order()
            {
                return mIGSpaceInterpolationOrder;
            }

//------------------------------------------------------------------------------

              virtual const Cluster * get_clusters_by_index( moris_index aCellClusterIndex ) const = 0;

//------------------------------------------------------------------------------

              virtual uint get_num_vertieces_on_set( const bool aOnlyPrimary ) = 0;

//------------------------------------------------------------------------------

              virtual moris::Matrix< DDSMat > get_vertieces_inds_on_block( const bool aOnlyPrimary ) = 0;

//------------------------------------------------------------------------------

              virtual uint get_num_cells_on_set( const bool aOnlyPrimary ) = 0;

//------------------------------------------------------------------------------

              virtual moris::Matrix< DDSMat > get_cell_inds_on_block( const bool aOnlyPrimary ) = 0;

//------------------------------------------------------------------------------
              virtual moris::uint
              get_num_clusters_on_set() const = 0;

//------------------------------------------------------------------------------

              virtual moris::Cell<Cluster const *>
              get_clusters_on_set() const = 0;

//------------------------------------------------------------------------------

              mtk::Geometry_Type get_auto_side_geometry_type( const mtk::Geometry_Type aGeometryType )
              {
                  mtk::Geometry_Type tSideGeometryType;

                  // depending on the parent geometry type
                  switch ( aGeometryType )
                  {
                      case ( mtk::Geometry_Type::LINE ):
                      {
                          tSideGeometryType = mtk::Geometry_Type::POINT;
                          break;
                      }
                      case ( mtk::Geometry_Type::QUAD ):
                      {
                          tSideGeometryType = mtk::Geometry_Type::LINE;
                          break;
                      }
                      case ( mtk::Geometry_Type::HEX ):
                      {
                          tSideGeometryType = mtk::Geometry_Type::QUAD;
                          break;
                      }
                      case ( mtk::Geometry_Type::TRI ):
                          {
                              tSideGeometryType = mtk::Geometry_Type::LINE;
                              break;
                          }
                      case ( mtk::Geometry_Type::TET ):
                          {
                              tSideGeometryType = mtk::Geometry_Type::TRI;
                              break;
                          }
                      default:
                      {
//                          MORIS_ERROR( false, " Geometry_Interpolator::get_auto_side_geometry_type - undefined geometry type. " );
                          tSideGeometryType = mtk::Geometry_Type::UNDEFINED;
                      }
                  }
                  return tSideGeometryType;
              }

//------------------------------------------------------------------------------

              mtk::Interpolation_Order get_auto_interpolation_order( const moris::uint        aNumVertices,
                                                                          const mtk::Geometry_Type aGeometryType )
              {
                  switch( aGeometryType )
                  {
                      case( mtk::Geometry_Type::LINE ) :
                          switch( aNumVertices )
                          {
                             case( 1 ) :
                                 return mtk::Interpolation_Order::UNDEFINED;
                                 break;
                             case( 2 ) :
                                 return mtk::Interpolation_Order::LINEAR;
                                 break;

                             case( 3 ) :
                                 return mtk::Interpolation_Order::QUADRATIC;
                                 break;

                             default :
                                 MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                                 return mtk::Interpolation_Order::UNDEFINED;
                                 break;
                          }

                      case( mtk::Geometry_Type::QUAD ) :
                          switch( aNumVertices )
                          {
                              case( 4 ) :
                                  return mtk::Interpolation_Order::LINEAR;
                                  break;

                              case( 8 ) :
                                  return mtk::Interpolation_Order::SERENDIPITY;
                                  break;

                              case( 9 ) :
                                  return mtk::Interpolation_Order::QUADRATIC;
                                  break;

                              case( 16 ) :
                                  return mtk::Interpolation_Order::CUBIC;
                                  break;

                              default :
                                  MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                                  return mtk::Interpolation_Order::UNDEFINED;
                                  break;
                          }

                      case( mtk::Geometry_Type::HEX ) :
                          switch( aNumVertices )
                          {
                              case( 8 ) :
                                  return mtk::Interpolation_Order::LINEAR;
                                  break;

                              case( 20 ) :
                                  return mtk::Interpolation_Order::SERENDIPITY;
                                  break;

                              case( 27 ) :
                                  return mtk::Interpolation_Order::QUADRATIC;
                                  break;

                              case( 64 ) :
                                  return mtk::Interpolation_Order::CUBIC;
                                  break;

                              default :
                                  MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                                  return mtk::Interpolation_Order::UNDEFINED;
                                  break;
                          }

                          case( mtk::Geometry_Type::TET ) :
                          switch( aNumVertices )
                          {
                              case( 4 ) :
                                  return mtk::Interpolation_Order::LINEAR;
                                  break;

                              case( 10 ) :
                                  return mtk::Interpolation_Order::QUADRATIC;
                                  break;

                              case( 20 ) :
                                  return mtk::Interpolation_Order::CUBIC;
                                  break;

                              default :
                                  MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                                  return mtk::Interpolation_Order::UNDEFINED;
                                  break;
                          }


                      default :
                          MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. ");
                          return mtk::Interpolation_Order::UNDEFINED;
                          break;
                  }
              }


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_SET_HPP_ */
