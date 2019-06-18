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

            // interpolation mesh geometry type
            mtk::Geometry_Type mIPGeometryType = mtk::Geometry_Type::UNDEFINED;

//------------------------------------------------------------------------------

        protected :

            moris::Cell<Cluster const *> mSetClusters;

            // integration mesh geometry type
            mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

//------------------------------------------------------------------------------

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

                MORIS_ASSERT( mIPGeometryType != mtk::Geometry_Type::UNDEFINED, " communicate_type(); undefined geometry type on all processors");
            };

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Set()
            { };

            Set(moris::Cell<Cluster const *>  aBlockSetClusters) : mSetClusters( aBlockSetClusters )
            {
                this->communicate_ip_geometry_type();
            };

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Set(){};

//------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */
//              virtual const moris::Matrix< DDUMat > &
//              get_list_of_block_cell_clusters() const = 0;

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

              virtual const Cluster  *
              get_cell_clusters_by_index( moris_index aCellClusterIndex ) const
              {
                  MORIS_ASSERT(false, "get_cell_clusters_by_index() virtual base class used");
                  return nullptr;
              };

              virtual const Cluster  *
              get_side_clusters_by_index( moris_index aCellClusterIndex ) const
              {
                  MORIS_ASSERT(false, "get_side_clusters_by_index() virtual base class used");
                  return nullptr;
              };

//------------------------------------------------------------------------------

              virtual const uint
              get_num_vertieces_on_set() const = 0;

//------------------------------------------------------------------------------

              virtual moris::Matrix< DDSMat >
              get_vertieces_inds_on_block() const = 0;

//------------------------------------------------------------------------------

              virtual const moris::uint
              get_num_clusters_on_set() const = 0;

//------------------------------------------------------------------------------

              virtual moris::Cell<Cluster const *>
              get_clusters_on_set() const
              {
                  MORIS_ASSERT(false, "get_cell_clusters_on_set() virtual base class used");
                  return moris::Cell<Cluster const *>(0);
              }

//------------------------------------------------------------------------------

              virtual moris::Cell<Cluster const *>
              get_side_clusters_on_set() const
              {
                  MORIS_ASSERT(false, "get_side_clusters_on_set() virtual base class used");
                  return moris::Cell<Cluster const *>(0);
              }

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
                          MORIS_ERROR( false, " Geometry_Interpolator::get_auto_side_geometry_type - undefined geometry type. " );
                          tSideGeometryType = mtk::Geometry_Type::UNDEFINED;
                      }
                  }
                  return tSideGeometryType;
              }


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_SET_HPP_ */
