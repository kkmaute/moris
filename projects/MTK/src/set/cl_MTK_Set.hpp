/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Set.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_SET_HPP_
#define SRC_MESH_CL_MTK_SET_HPP_

#include <string>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "fn_unique.hpp"         //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK/src
#include "cl_MTK_Cell.hpp"      //MTK/src

#include "cl_MTK_Cell_Cluster.hpp"       //MTK/src
#include "cl_MTK_Side_Cluster.hpp"       //MTK/src
#include "cl_Communication_Tools.hpp"    //MTK/src

#include "cl_MTK_Set_Communicator.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        class Set
        {
            friend mtk::Set_Communicator;

            //------------------------------------------------------------------------------

          private:
            // name of the set
            std::string mSetName;

            // space interpolation order for IP cells
            mtk::Interpolation_Order mIPSpaceInterpolationOrder;

            // space interpolation order for IG cells
            mtk::Interpolation_Order mIGSpaceInterpolationOrder;

            moris_index mSetIndex = MORIS_INDEX_MAX;

            mtk::CellShape mIGCellShape = CellShape::UNDEFINED;
            mtk::CellShape mIPCellShape = CellShape::UNDEFINED;

            bool mIsTrivialLeader   = false;
            bool mIsTrivialFollower = false;
            bool mLeaderLock        = false;
            bool mFollowerLock      = false;

            mtk::CellTopology mCellTopology = CellTopology::UNDEFINED;

            moris::uint mSpatialDim;

            Matrix< IndexMat > mSetColors;

            //------------------------------------------------------------------------------

          protected:
            Vector< Cluster const * > mSetClusters;

            mtk::Geometry_Type mIPGeometryType = mtk::Geometry_Type::UNDEFINED;
            mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

            SetType mSetType = SetType::UNDEFINED;

            bool mOwnedByPeriodicBCFlag = false;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Set(){};

            //------------------------------------------------------------------------------

            Set(
                    std::string const               &aName,
                    Vector< Cluster const * > const &aBlockSetClusters,
                    Matrix< IndexMat > const        &aColors,
                    uint const                      &aSpatialDim )
                    : mSetName( aName )
                    , mSpatialDim( aSpatialDim )
                    , mSetColors( aColors )
                    , mSetClusters( aBlockSetClusters )
            {
                this->init_ip_geometry_type();
                this->init_interpolation_order();
            }

            //------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual ~Set(){};

            //------------------------------------------------------------------------------

            std::string
            get_set_name() const
            {
                return mSetName;
            }

            //------------------------------------------------------------------------------

            SetType
            get_set_type() const
            {
                return mSetType;
            }

            //------------------------------------------------------------------------------

            uint
            get_clusters() const
            {
                MORIS_ERROR( !( mSpatialDim < 1 ) || !( mSpatialDim > 3 ), "Set::get_spatial_dim(), Spatial dim < 1 or > 3" );
                return mSpatialDim;
            }

            //------------------------------------------------------------------------------

            uint
            get_spatial_dim() const
            {
                MORIS_ERROR( !( mSpatialDim < 1 ) || !( mSpatialDim > 3 ), "Set::get_spatial_dim(), Spatial dim < 1 or > 3" );
                return mSpatialDim;
            }

            //------------------------------------------------------------------------------

            void
            set_set_index( moris_index aIndex )
            {
                mSetIndex = aIndex;
            }

            //------------------------------------------------------------------------------

            moris_index
            get_set_index() const
            {
                MORIS_ASSERT( mSetIndex != MORIS_INDEX_MAX, "Set::get_set_index(), Set index not set" );
                return mSetIndex;
            }

            //------------------------------------------------------------------------------

            void
            set_cell_topology( CellTopology aCellTopology )
            {
                mCellTopology = aCellTopology;
            }

            //------------------------------------------------------------------------------

            CellTopology
            get_cell_topology() const
            {
                MORIS_ASSERT( mCellTopology != CellTopology::UNDEFINED, "Set::get_cell_topology(), Cell topology not set" );
                return mCellTopology;
            }

            //------------------------------------------------------------------------------

            void
            set_IG_cell_shape( CellShape aCellShape )
            {
                mIGCellShape = aCellShape;
            }

            //------------------------------------------------------------------------------

            void
            set_IP_cell_shape( CellShape aCellShape )
            {
                mIPCellShape = aCellShape;
            }

            //------------------------------------------------------------------------------

            CellShape
            get_IG_cell_shape() const
            {
                MORIS_ASSERT( mIGCellShape != CellShape::UNDEFINED, "Set::get_IG_cell_shape(), Cell shape not set" );
                return mIGCellShape;
            }

            //------------------------------------------------------------------------------

            CellShape
            get_IP_cell_shape() const
            {
                MORIS_ASSERT( mIPCellShape != CellShape::UNDEFINED, "Set::get_IG_cell_shape(), Cell shape not set" );
                return mIPCellShape;
            }

            //------------------------------------------------------------------------------

            Matrix< IndexMat > const &
            get_set_colors() const
            {
                return mSetColors;
            }
            //------------------------------------------------------------------------------

            /**
             * return a label that describes the block
             */

            bool
            is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER )
            {
                if ( !mLeaderLock && aIsLeader == mtk::Leader_Follower::LEADER )
                {
                    this->communicate_is_trivial_flag( aIsLeader );

                    mLeaderLock = true;
                }
                else if ( !mFollowerLock && aIsLeader == mtk::Leader_Follower::FOLLOWER )
                {
                    this->communicate_is_trivial_flag( aIsLeader );

                    mFollowerLock = true;
                }

                if ( aIsLeader == mtk::Leader_Follower::LEADER )
                {
                    return mIsTrivialLeader;
                }
                else if ( aIsLeader == mtk::Leader_Follower::FOLLOWER )
                {
                    return mIsTrivialFollower;
                }
                else
                {
                    MORIS_ASSERT( false, " is_trivial(); undefined type. has to be leader or follower" );

                    return false;
                }
            }

            //------------------------------------------------------------------------------
            //------------------------------------------------------------------------------

            mtk::Geometry_Type
            get_interpolation_cell_geometry_type() const
            {
                return mIPGeometryType;
            }

            //------------------------------------------------------------------------------

            mtk::Geometry_Type
            get_integration_cell_geometry_type() const
            {
                return mIGGeometryType;
            }

            //------------------------------------------------------------------------------

            mtk::Interpolation_Order
            get_interpolation_cell_interpolation_order() const
            {
                return mIPSpaceInterpolationOrder;
            }

            //------------------------------------------------------------------------------

            mtk::Interpolation_Order
            get_integration_cell_interpolation_order() const
            {
                return mIGSpaceInterpolationOrder;
            }

            //------------------------------------------------------------------------------

            virtual const Cluster *get_clusters_by_index( moris_index aCellClusterIndex ) const = 0;

            //------------------------------------------------------------------------------

            virtual uint get_num_vertices_on_set( bool aOnlyPrimary ) = 0;

            //------------------------------------------------------------------------------

            virtual moris::Matrix< DDSMat > get_ig_vertices_inds_on_block( bool aOnlyPrimary ) = 0;

            //------------------------------------------------------------------------------

            virtual uint get_num_cells_on_set( bool aOnlyPrimary ) = 0;

            //------------------------------------------------------------------------------

            virtual moris::Matrix< DDSMat > get_cell_inds_on_block( bool aOnlyPrimary ) = 0;

            //------------------------------------------------------------------------------
            virtual moris::uint
            get_num_clusters_on_set() const = 0;

            //------------------------------------------------------------------------------

            virtual Vector< Cluster const * >
            get_clusters_on_set() const = 0;

            //------------------------------------------------------------------------------

            mtk::Geometry_Type
            get_auto_side_geometry_type( const mtk::Geometry_Type aGeometryType ) const
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

            mtk::Interpolation_Order
            get_auto_interpolation_order(
                    const moris::uint        aNumVertices,
                    const mtk::Geometry_Type aGeometryType ) const
            {
                switch ( aGeometryType )
                {
                    case mtk::Geometry_Type::LINE:
                    {
                        switch ( aNumVertices )
                        {
                            case 1:
                                return mtk::Interpolation_Order::UNDEFINED;
                            case 2:
                                return mtk::Interpolation_Order::LINEAR;

                            case 3:
                                return mtk::Interpolation_Order::QUADRATIC;

                            default:
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. " );
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::QUAD:
                    {
                        switch ( aNumVertices )
                        {
                            case 4:
                                return mtk::Interpolation_Order::LINEAR;

                            case 8:
                                return mtk::Interpolation_Order::SERENDIPITY;

                            case 9:
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 16:
                                return mtk::Interpolation_Order::CUBIC;

                            default:
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. " );
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::HEX:
                    {
                        switch ( aNumVertices )
                        {
                            case 8:
                                return mtk::Interpolation_Order::LINEAR;

                            case 20:
                                return mtk::Interpolation_Order::SERENDIPITY;

                            case 27:
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 64:
                                return mtk::Interpolation_Order::CUBIC;

                            default:
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. " );
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }
                    case mtk::Geometry_Type::TET:
                    {
                        switch ( aNumVertices )
                        {
                            case 4:
                                return mtk::Interpolation_Order::LINEAR;

                            case 10:
                                return mtk::Interpolation_Order::QUADRATIC;

                            case 20:
                                return mtk::Interpolation_Order::CUBIC;

                            default:
                                MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TET and number of vertices. " );
                                return mtk::Interpolation_Order::UNDEFINED;
                        }
                    }

                    default:
                        MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. " );
                        return mtk::Interpolation_Order::UNDEFINED;
                }
            }

            //-----------------------------------------------------------------------------

            bool
            give_ownership_periodic()
            {
                mOwnedByPeriodicBCFlag = true;

                return mOwnedByPeriodicBCFlag;
            }

            //-----------------------------------------------------------------------------

            virtual size_t
            capacity()
            {
                size_t tTotalSize = 0;

                // name of the set
                tTotalSize += sizeof( mSetName );
                tTotalSize += sizeof( mIPGeometryType );
                tTotalSize += sizeof( mIPSpaceInterpolationOrder );
                tTotalSize += sizeof( mIGSpaceInterpolationOrder );
                tTotalSize += sizeof( mIsTrivialLeader );
                tTotalSize += sizeof( mIsTrivialFollower );
                tTotalSize += sizeof( mLeaderLock );
                tTotalSize += sizeof( mFollowerLock );
                tTotalSize += sizeof( mSetIndex );
                tTotalSize += sizeof( mIGCellShape );
                tTotalSize += sizeof( mCellTopology );
                tTotalSize += sizeof( mIPCellShape );
                tTotalSize += sizeof( mSpatialDim );
                tTotalSize += sizeof( mSetColors ) + mSetColors.capacity();
                tTotalSize += sizeof( mIGGeometryType );
                tTotalSize += sizeof( mSetType );
                tTotalSize += sizeof( mOwnedByPeriodicBCFlag );
                tTotalSize += mSetClusters.capacity() * ( 1 + sizeof( void * ) );

                return tTotalSize;
            }

            //-----------------------------------------------------------------------------

          protected:
            //-----------------------------------------------------------------------------

            void
            init_ip_geometry_type()
            {
                mIPGeometryType = mtk::Geometry_Type::UNDEFINED;

                if ( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    mIPGeometryType = mSetClusters( 0 )->get_interpolation_cell().get_geometry_type();
                }

                // TODO: check if it works with this commented out
                // uint tRecIPGeometryType = min_all( (uint)mIPGeometryType );
                // mIPGeometryType = static_cast< mtk::Geometry_Type >( tRecIPGeometryType );

                // MORIS_ASSERT( mIPGeometryType != mtk::Geometry_Type::UNDEFINED, " communicate_type(); undefined geometry type on all processors");
            }

            //------------------------------------------------------------------------------

            // FIXME should be user-defined in FEM
            void
            init_interpolation_order()
            {
                // MORIS_ASSERT(
                //         mIPGeometryType != mtk::Geometry_Type::UNDEFINED,
                //         " init_interpolation_order(); undefined geometry type on this processor. Try calling communicate_ip_geometry_type() first.");

                mIPSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;
                mIGSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

                if ( mSetClusters.size() > 0 )
                {
                    // interpolation order for IP cells fixme
                    mIPSpaceInterpolationOrder = mSetClusters( 0 )->get_interpolation_cell( mtk::Leader_Follower::LEADER ).get_interpolation_order();

                    // get list of primary IG cells in cluster
                    Vector< mtk::Cell const * > const &tPrimaryIgCellsInCluster = mSetClusters( 0 )->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );

                    // set interpolation order for IG cells fixme
                    if ( tPrimaryIgCellsInCluster.size() > 0 )
                    {
                        mIGSpaceInterpolationOrder = tPrimaryIgCellsInCluster( 0 )->get_interpolation_order();
                    }
                    else    // in case there are void clusters, look at the void clusters
                    {
                        mIGSpaceInterpolationOrder = mSetClusters( 0 )->get_void_cells_in_cluster()( 0 )->get_interpolation_order();
                    }
                }

                // TODO: check if it works with this commented out
                // uint tRecIPInterpolationOrder = min_all( (uint)mIPSpaceInterpolationOrder );
                // uint tRecIGInterpolationOrder = min_all( (uint)mIGSpaceInterpolationOrder );
                // mIPSpaceInterpolationOrder = static_cast< mtk::Interpolation_Order >( tRecIPInterpolationOrder );
                // mIGSpaceInterpolationOrder = static_cast< mtk::Interpolation_Order >( tRecIGInterpolationOrder );

                // MORIS_ASSERT( mIPSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " init_interpolation_order(); undefined ip interpolation order on this processor");
                // MORIS_ASSERT( mIGSpaceInterpolationOrder != mtk::Interpolation_Order::UNDEFINED, " init_interpolation_order(); undefined ig interpolation order on this processor");
            }

            //------------------------------------------------------------------------------

            // FIXME should be user-defined in FEM
            void
            communicate_is_trivial_flag( const mtk::Leader_Follower aIsLeader )
            {
                sint tIsTrivial = 1;

                if ( mSetClusters.size() > 0 )
                {
                    // set the integration geometry type
                    tIsTrivial = (sint)mSetClusters( 0 )->is_trivial( aIsLeader );    // FIXME change for double sided set
                }

                sint tIsTrivialMax = max_all( tIsTrivial );

                if ( tIsTrivialMax == 1 && aIsLeader == mtk::Leader_Follower::LEADER )
                {
                    mIsTrivialLeader = true;
                }
                else if ( tIsTrivialMax == 1 && aIsLeader == mtk::Leader_Follower::FOLLOWER )
                {
                    mIsTrivialFollower = true;
                }
            }

            //------------------------------------------------------------------------------
            //------------------------------------------------------------------------------

            void
            set_interpolation_cell_geometry_type( mtk::Geometry_Type aGeometryType )
            {
                MORIS_ASSERT(
                        mIPGeometryType == mtk::Geometry_Type::UNDEFINED || mIPGeometryType == aGeometryType,
                        "mtk::Set::set_interpolation_cell_geometry_type() - "
                        "IP cell geometry type is already set and is different. This is a conflict." );
                mIPGeometryType = aGeometryType;
            }

            //------------------------------------------------------------------------------

            void
            set_integration_cell_geometry_type( mtk::Geometry_Type aGeometryType )
            {
                MORIS_ASSERT(
                        mIGGeometryType == mtk::Geometry_Type::UNDEFINED || mIGGeometryType == aGeometryType,
                        "mtk::Set::set_integration_cell_geometry_type() - "
                        "IG cell geometry type is already set and is different. This is a conflict." );
                mIGGeometryType = aGeometryType;
            }

            //------------------------------------------------------------------------------

            void
            set_interpolation_cell_interpolation_order( mtk::Interpolation_Order aIpOrder )
            {
                MORIS_ASSERT(
                        mIPSpaceInterpolationOrder == mtk::Interpolation_Order::UNDEFINED || mIPSpaceInterpolationOrder == aIpOrder,
                        "mtk::Set::set_interpolation_cell_interpolation_order() - "
                        "IP cell interpolation order is already set and is different. This is a conflict." );
                mIPSpaceInterpolationOrder = aIpOrder;
            }

            //------------------------------------------------------------------------------

            void
            set_integration_cell_interpolation_order( mtk::Interpolation_Order aIpOrder )
            {
                MORIS_ASSERT(
                        mIGSpaceInterpolationOrder == mtk::Interpolation_Order::UNDEFINED || mIGSpaceInterpolationOrder == aIpOrder,
                        "mtk::Set::set_integration_cell_interpolation_order() - "
                        "IG cell interpolation order is already set and is different. This is a conflict." );
                mIGSpaceInterpolationOrder = aIpOrder;
            }

            //------------------------------------------------------------------------------

        };    // end class: mtk::Set

        //------------------------------------------------------------------------------

    } /* namespace mtk */
} /* namespace moris */

//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_SET_HPP_ */
