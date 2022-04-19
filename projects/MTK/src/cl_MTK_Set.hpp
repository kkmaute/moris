/*
 * cl_MTK_Element_Set.hpp
 *
 *  Created on: Jun 24, 2019
 *      Author: schmidt
 */

#ifndef SRC_MESH_CL_MTK_SET_HPP_
#define SRC_MESH_CL_MTK_SET_HPP_

#include <string>

#include "typedefs.hpp"       //MRS/COR/src
#include "fn_unique.hpp"      //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp"  //MTK/src
#include "cl_MTK_Cell.hpp"    //MTK/src

#include "cl_MTK_Cell_Cluster.hpp"     //MTK/src
#include "cl_MTK_Side_Cluster.hpp"     //MTK/src
#include "cl_Communication_Tools.hpp"  //MTK/src

#include "cl_Mesh_Enums.hpp"

namespace moris
{
    namespace mtk
    {
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

            moris_index mSetIndex = MORIS_INDEX_MAX;

            enum CellTopology mCellTopology = CellTopology::END_ENUM;

            enum CellShape mIGCellShape = CellShape::END_ENUM;

            enum CellShape mIPCellShape = CellShape::END_ENUM;

            moris::uint mSpatialDim;

            Matrix<IndexMat> mSetColors;


            //------------------------------------------------------------------------------

        protected :

            moris::Cell< Cluster const * > mSetClusters;

            // integration mesh geometry type
            mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

            moris::SetType mSetType = moris::SetType::END_ENUM;

            bool                          mOwendbyPeriodicBCFlag = false;

        public:

            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Set() { };

            //------------------------------------------------------------------------------

            Set(
                    std::string                  const & aName,
                    moris::Cell<Cluster const *> const & aBlockSetClusters,
                    Matrix<IndexMat>             const & aColors,
                    uint                         const & aSpatialDim )
                    : mSetName( aName ),
                      mSpatialDim( aSpatialDim ),
                      mSetColors(aColors),
                      mSetClusters( aBlockSetClusters )
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

            enum moris::SetType get_set_type()
            {
                return mSetType;
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

            void set_set_index( moris_index aIndex )
            {
                mSetIndex = aIndex;
            }

            //------------------------------------------------------------------------------

            moris_index get_set_index(  )
            {
                MORIS_ASSERT( mSetIndex != MORIS_INDEX_MAX, "Set::get_set_index(), Set index not set" );
                return mSetIndex;
            }

            //------------------------------------------------------------------------------

            void set_cell_topology( enum CellTopology aCellTopology )
            {
                mCellTopology = aCellTopology;
            }

            //------------------------------------------------------------------------------

            enum CellTopology get_cell_topology()
            {
                MORIS_ASSERT( mCellTopology != CellTopology::END_ENUM, "Set::get_cell_topology(), Cell topology not set" );
                return mCellTopology;
            }

            
            //------------------------------------------------------------------------------

            void set_IG_cell_shape( enum CellShape aCellShape )
            {
                mIGCellShape = aCellShape;
            }

            //------------------------------------------------------------------------------

            void set_IP_cell_shape( enum CellShape aCellShape )
            {
                mIPCellShape = aCellShape;
            }

            //------------------------------------------------------------------------------

            enum CellShape get_IG_cell_shape()
            {
                MORIS_ASSERT( mIGCellShape != CellShape::END_ENUM, "Set::get_IG_cell_shape(), Cell shape not set" );
                return mIGCellShape;
            }

            //------------------------------------------------------------------------------

            enum CellShape get_IP_cell_shape()
            {
                MORIS_ASSERT( mIPCellShape != CellShape::END_ENUM, "Set::get_IG_cell_shape(), Cell shape not set" );
                return mIPCellShape;
            }

            //------------------------------------------------------------------------------

            Matrix<IndexMat> const &
            get_set_colors()
            {
                return mSetColors;
            }
            //------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */

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

            virtual uint get_num_vertices_on_set( const bool aOnlyPrimary ) = 0;

            //------------------------------------------------------------------------------

            virtual moris::Matrix< DDSMat > get_ig_vertices_inds_on_block( const bool aOnlyPrimary ) = 0;

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
                    case mtk::Geometry_Type::LINE:
                    {
                        tSideGeometryType = mtk::Geometry_Type::POINT;
                        break;
                    }
                    case mtk::Geometry_Type::QUAD:
                    {
                        tSideGeometryType = mtk::Geometry_Type::LINE;
                        break;
                    }
                    case mtk::Geometry_Type::HEX:
                    {
                        tSideGeometryType = mtk::Geometry_Type::QUAD;
                        break;
                    }
                    case mtk::Geometry_Type::TRI:
                    {
                        tSideGeometryType = mtk::Geometry_Type::LINE;
                        break;
                    }
                    case mtk::Geometry_Type::TET:
                    {
                        tSideGeometryType = mtk::Geometry_Type::TRI;
                        break;
                    }
                    default:
                    {
                        // MORIS_ERROR( false, " Geometry_Interpolator::get_auto_side_geometry_type - undefined geometry type. " );
                        tSideGeometryType = mtk::Geometry_Type::UNDEFINED;
                    }
                }
                return tSideGeometryType;
            }

            //------------------------------------------------------------------------------

            mtk::Interpolation_Order get_auto_interpolation_order(
                    const moris::uint        aNumVertices,
                    const mtk::Geometry_Type aGeometryType )
            {
                switch( aGeometryType )
                {
                    case mtk::Geometry_Type::LINE:
                    {
                        switch( aNumVertices )
                        {
                            case 1 :
                                return mtk::Interpolation_Order::UNDEFINED;
                            case 2 :
                                return mtk::Interpolation_Order::LINEAR;

                            case 3 :
                                return mtk::Interpolation_Order::QUADRATIC;

                            default :
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::QUAD:
                    {
                        switch( aNumVertices )
                        {
                            case 4 :
                                return mtk::Interpolation_Order::LINEAR;

                            case 8 :
                                return mtk::Interpolation_Order::SERENDIPITY;

                            case 9 :
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 16 :
                                return mtk::Interpolation_Order::CUBIC;

                            default :
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::HEX :
                    {
                        switch( aNumVertices )
                        {
                            case 8 :
                                return mtk::Interpolation_Order::LINEAR;

                            case 20 :
                                return mtk::Interpolation_Order::SERENDIPITY;

                            case 27 :
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 64 :
                                return mtk::Interpolation_Order::CUBIC;

                            default :
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::TET :
                    {
                        switch( aNumVertices )
                        {
                            case 4 :
                                return mtk::Interpolation_Order::LINEAR;

                            case 10 :
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 20 :
                                return mtk::Interpolation_Order::CUBIC;

                            default :
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TET and number of vertices. ");
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }

                    default :
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. ");
                        return mtk::Interpolation_Order::UNDEFINED;
                }
            }

            //-----------------------------------------------------------------------------

            bool
            give_ownership_periodic()
            {
                 mOwendbyPeriodicBCFlag = true;

                return  mOwendbyPeriodicBCFlag;
            }

            //-----------------------------------------------------------------------------

            size_t
            capacity()
            {
                size_t tTotalSize = 0;

                // name of the set
                tTotalSize += sizeof( mSetName );
                tTotalSize += sizeof( mIPGeometryType );
                tTotalSize += sizeof( mIPSpaceInterpolationOrder );
                tTotalSize += sizeof( mIGSpaceInterpolationOrder );
                tTotalSize += sizeof( mIsTrivialMaster );
                tTotalSize += sizeof( mIsTrivialSlave );
                tTotalSize += sizeof( mMasterLock );
                tTotalSize += sizeof( mSlaveLock );
                tTotalSize += sizeof( mSetIndex );
                tTotalSize += sizeof( mIGCellShape );
                tTotalSize += sizeof( mCellTopology );
                tTotalSize += sizeof( mIPCellShape );
                tTotalSize += sizeof( mSpatialDim );
                tTotalSize += sizeof( mSetColors ) + mSetColors.capacity();
                tTotalSize += sizeof( mIGGeometryType );
                tTotalSize += sizeof( mSetType );
                tTotalSize += sizeof( mOwendbyPeriodicBCFlag );
                tTotalSize += mSetClusters.capacity() * ( 1 + sizeof( void * ) );

                return tTotalSize;
            }

            //-----------------------------------------------------------------------------

        protected:

            void communicate_ip_geometry_type()
            {
                mtk::Geometry_Type tIPGeometryType = mtk::Geometry_Type::UNDEFINED;

                if( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    tIPGeometryType = mSetClusters( 0 )->get_interpolation_cell().get_geometry_type();
                }

                uint tRecIPGeometryType = min_all( (uint)tIPGeometryType );

                mIPGeometryType = static_cast<enum mtk::Geometry_Type> (tRecIPGeometryType);

                //                MORIS_ASSERT( mIPGeometryType != mtk::Geometry_Type::UNDEFINED, " communicate_type(); undefined geometry type on all processors");
            };

            //------------------------------------------------------------------------------

            // FIXME should be user-defined in FEM
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

                    // get list of primary IG cells in cluster
                    moris::Cell< mtk::Cell const* > const& tPrimaryIgCellsInCluster = mSetClusters( 0 )->get_primary_cells_in_cluster( mtk::Master_Slave::MASTER );

                    // set interpolation order for IG cells fixme
                    if( tPrimaryIgCellsInCluster.size() > 0 )
                    {
                        tIGInterpolationOrder = tPrimaryIgCellsInCluster( 0 )->get_interpolation_order();
                    }
                    else // in case there are void clusters, look at the void clusters
                    {
                        tIGInterpolationOrder = mSetClusters( 0 )->get_void_cells_in_cluster()( 0 )->get_interpolation_order();
                    }
                }

                uint tRecIPInterpolationOrder = min_all( (uint)tIPInterpolationOrder );
                uint tRecIGInterpolationOrder = min_all( (uint)tIGInterpolationOrder );

                mIPSpaceInterpolationOrder = static_cast<enum mtk::Interpolation_Order> (tRecIPInterpolationOrder);
                mIGSpaceInterpolationOrder = static_cast<enum mtk::Interpolation_Order> (tRecIGInterpolationOrder);

                //                MORIS_ASSERT( mIPSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " communicate_interpolation_order(); undefined ip interpolation order on this processor");
                //                MORIS_ASSERT( mIGSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " communicate_interpolation_order(); undefined ig interpolation order on this processor");
            }

            //------------------------------------------------------------------------------

            // FIXME should be user-defined in FEM
            void communicate_is_trivial_flag( const mtk::Master_Slave aIsMaster )
            {
                sint tIsTrivial = 1;

                if( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    tIsTrivial = (sint)mSetClusters( 0 )->is_trivial( aIsMaster ); //FIXME change for double sided set
                }

                sint tIsTrivialMax = max_all( tIsTrivial );

                if( tIsTrivialMax == 1 && aIsMaster == mtk::Master_Slave::MASTER )
                {
                    mIsTrivialMaster = true;
                }
                else if( tIsTrivialMax == 1 && aIsMaster == mtk::Master_Slave::SLAVE )
                {
                    mIsTrivialSlave = true;
                }
            };
        };
    } /* namespace mtk */
} /* namespace moris */

//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_SET_HPP_ */
