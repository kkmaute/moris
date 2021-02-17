/*
 * cl_FEM_Node.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: messe
 */

#ifndef PROJECTS_FEM_SRC_CL_FEM_NODE_HPP_
#define PROJECTS_FEM_SRC_CL_FEM_NODE_HPP_

#include "typedefs.hpp"           //MRS/COR/src
#include "cl_MTK_Vertex.hpp"      //MTK/src
#include "cl_FEM_Node_Base.hpp"   //MTK/src

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class Node : public Node_Base
        {
            private:

                // mtk vertex pointer
                const mtk::Vertex * mVertex;

                // storage for xyz pdv active flags
                Matrix< DDSMat > mIsActiveXYZ;

                // storage for xyz pdv ids
                Matrix< DDSMat > mXYZPdvIds;

                // storage for xyz pdv local cluster assembly indices
                Matrix< DDSMat > mXYZLocalAssemblyIndices;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                Node( const mtk::Vertex * aVertex ) : mVertex( aVertex )
                {
                    //                mOwner = mVertex->get_owner();
                }

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Node(){};

                //------------------------------------------------------------------------------
                /**
                 * returns the owner of this node
                 */
                auto get_owner() const -> decltype( mVertex->get_owner() )
                {
                    return mVertex->get_owner();
                }

                //------------------------------------------------------------------------------

                /**
                 * returns the T-Matrix of this node
                 */
                const Matrix< DDRMat > * get_t_matrix( const uint aBSplineMeshIndex ) const
                {
                    return mVertex->get_interpolation( aBSplineMeshIndex )->get_weights();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the B-Spline IDs of this node
                 */
                Matrix< IdMat > get_adof_ids( const uint aBSplineMeshIndex ) const
                {
                    return mVertex->get_interpolation( aBSplineMeshIndex )->get_ids();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the B-Spline IDs of this node
                 */
                Matrix< IndexMat > get_adof_indices( const uint aBSplineMeshIndex ) const
                {
                    return mVertex->get_interpolation( aBSplineMeshIndex )->get_indices();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the proc owners of the IDs of this node
                 */
                Matrix< IdMat > get_adof_owners( const uint aBSplineMeshIndex ) const
                {
                    return mVertex->get_interpolation( aBSplineMeshIndex )->get_owners();
                }

                //------------------------------------------------------------------------------
                /**
                 * get the ID of this node
                 * @param[ in ] id id for this node
                 */
                moris_id get_id() const
                {
                    return mVertex->get_id();
                }

                //------------------------------------------------------------------------------
                /**
                 * get the Index of this node
                 * @param[ out ] index index for this node
                 */
                moris_index get_index() const
                {
                    return mVertex->get_index();
                }

                //------------------------------------------------------------------------------
                /**
                 * is owned
                 * @param[ out ] bool true if node is owned by processor
                 */
                bool id_owned()
                {
                    bool tOwned = true;

                    if( mVertex->get_owner() != par_rank() )
                    {
                        tOwned = false;
                    }

                    return tOwned;
                }

                //------------------------------------------------------------------------------
                /**
                 * get vertex coordinates
                 * @param[ in ] aVertexCoords matrix to fill with vertex coordinates
                 */
                void get_vertex_coords( Matrix< DDRMat > & aVertexCoords )
                {
                    aVertexCoords = mVertex->get_coords();
                }

                //------------------------------------------------------------------------------
                /**
                 * get vertex active flags (if relevant)
                 * @param[ out ] aPdvTypes
                 * @param[ out ] aIsActiveDv
                 */
                void get_vertex_xyz_active_flags(
                         Matrix< DDSMat >                   & aIsActiveDv,
                        const moris::Cell < enum PDV_Type > & aPdvTypes )
                {
                    // get number of requested pdv types
                    uint tNumPdvTypes = aPdvTypes.size();

                    // set size for active flag
                    aIsActiveDv.set_size( 1, tNumPdvTypes );

                    // loop over requested pdv types
                    for( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
                    {
                        // get pdv index
                        uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                        // set value
                        aIsActiveDv( iPdvType ) = mIsActiveXYZ( tXYZIndex );
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * set vertex active flags (if relevant)
                 * @param[ out ] aIsActiveDv
                 */
                void set_vertex_xyz_active_flags(
                        moris::Cell< Matrix< DDSMat > > & aIsActiveDv )
                {
                    // get num of pdv
                    uint tNumXYZPdv = aIsActiveDv.size();

                    // init mIsActiveXYZ
                    mIsActiveXYZ.set_size( 1, tNumXYZPdv, 0.0 );

                    // fill mIsActiveXYZ
                    for( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
                    {
                        mIsActiveXYZ( iPdvType ) = aIsActiveDv( iPdvType )( 0 );
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * set vertex active flags (if relevant)
                 * @param[ in ] aXYZPvIds list of xyz pdv ids
                 */
                void set_vertex_xyz_pdv_ids(
                        moris::Cell< Matrix< DDSMat > > & aXYZPvIds )
                {
                    // get num of pdv
                    uint tNumXYZPdv = aXYZPvIds.size();

                    // init mXYZPdvIds
                    mXYZPdvIds.set_size( 1, tNumXYZPdv, -1 );

                    // fill mXYZPdvIds
                    for( uint iPdvType = 0; iPdvType < tNumXYZPdv; iPdvType++ )
                    {
                        mXYZPdvIds( iPdvType ) = aXYZPvIds( iPdvType )( 0 );
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * get vertex xyz pdv ids (if relevant)
                 * @param[ in ] aXYZPdvIds
                 * @param[ in ] aIsActiveDv
                 */
                void get_vertex_xyz_pdv_ids(
                        Matrix< DDSMat >                    & aXYZPdvIds,
                        const moris::Cell < enum PDV_Type > & aPdvTypes )
                {
                    // get number of requested pdv types
                    uint tNumPdvTypes = aPdvTypes.size();

                    // set size for pdv ids
                    aXYZPdvIds.set_size( 1, tNumPdvTypes, -1 );

                    // loop over requested pdv types
                    for( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
                    {
                        // get pdv index
                        uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                        // set value
                        aXYZPdvIds( iPdvType ) = mXYZPdvIds( tXYZIndex );
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * reset xyz pdv local cluster assembly indices
                 */
                void reset_local_xyz_pdv_assembly_index()
                {
                    // get spatial dimension
                    uint tSpaceDim = mVertex->get_coords().numel();

                    // reset assembly indices on node
                    mXYZLocalAssemblyIndices.set_size( tSpaceDim, 1, -1 );
                }

                //------------------------------------------------------------------------------
                /**
                 * set x/y/z pdv local cluster assembly index
                 * @param[ in ] aLocalPdvAssemblyIndex a local cluster assembly index
                 * @param[ in ] aPdvType               enum for x/y/z pdv
                 */
                void set_local_xyz_pdv_assembly_index(
                        moris_index   aLocalPdvAssemblyIndex,
                        enum PDV_Type aPdvType )
                {
                    mXYZLocalAssemblyIndices( static_cast< uint >( aPdvType ) ) = aLocalPdvAssemblyIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * get x/y/z pdv local cluster assembly indices
                 * @param[ in ] aXYZLocalAssemblyIndices matrix to fill with local cluster assembly indices
                 * @param[ in ] aPdvType               list of enums for requested x/y/z pdv
                 */
                void get_local_xyz_pdv_assembly_indices(
                        Matrix< DDSMat >                    & aXYZLocalAssemblyIndices,
                        const moris::Cell < enum PDV_Type > & aPdvTypes )
                {
                     // get number of requested pdv types
                     uint tNumPdvTypes = aPdvTypes.size();

                     // set size for pdv ids
                     aXYZLocalAssemblyIndices.set_size( 1, tNumPdvTypes, -1 );

                     // loop over requested pdv types
                     for( uint iPdvType = 0; iPdvType < tNumPdvTypes; iPdvType++ )
                     {
                         // get pdv index
                         uint tXYZIndex = static_cast< uint >( aPdvTypes( iPdvType ) );

                         // set value
                         aXYZLocalAssemblyIndices( iPdvType ) = mXYZLocalAssemblyIndices( tXYZIndex );
                     }
                }

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace MSI */
} /* namespace moris */



#endif /* PROJECTS_FEM_SRC_CL_NODE_HPP_ */
