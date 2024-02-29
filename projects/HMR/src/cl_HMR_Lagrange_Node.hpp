/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Node.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_NODE_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_NODE_HPP_
#include "cl_HMR_Basis.hpp"
#include "cl_HMR_Lagrange_Node_Interpolation.hpp"
#include "moris_typedefs.hpp" //COR/src
#include "cl_HMR_Lagrange_Node_Interpolation.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------
    /**
     * \brief Lagrange Node class, templated against dimension
     */
    template< uint N >
    class Lagrange_Node : public Basis
    {
        //! local ijk position on proc
        luint         mIJK[ N ];

        //! global coordinates
        real          mXYZ[ N ];

        //! the T-Matrix of this node
        // Matrix< DDRMat >   mTMatrix;

        //! interpolator object
        Vector< Lagrange_Node_Interpolation * > mInterpolations;

        //! bitset telling if interpolation is set
        Bitset< gNumberOfMeshes > mHaveInterpolation;

        bool mHaveInterpolationContainer = false;

// ----------------------------------------------------------------------------
        public:
// ----------------------------------------------------------------------------

        /**
         * default constructor
         *
         * @param[in]   aIJK          ijk position of node
         * @param[in]   aLevel        level on which basis exists
         * @param[in]   aOwner        owner of basis
         *
         */
        Lagrange_Node(
                const luint        * aIJK,
                const  uint        & aLevel,
                const  uint        & aOwner) :
                    Basis( aLevel, aOwner ),
                    mInterpolations( gNumberOfMeshes, nullptr )
        {
            // save ijk position in memory.
            for( uint k=0; k<N; ++k )
            {
                mIJK[ k ] = aIJK[ k ];
            }
        }

// ----------------------------------------------------------------------------

        /**
         * default destructor
         */
        ~Lagrange_Node()
        {
            // delete element container
            if ( mNumberOfConnectedElements != 0 )
            {
                mElements.clear();
            }

            // delete facet container
            this->delete_facet_container();

            // delete edge container
            if( N == 3 )
            {
                this->delete_edge_container();
            }

            if( mHaveInterpolationContainer )
            {
                // delete interpolation objects
                for( uint k=0; k<gNumberOfMeshes; ++k )
                {
                    if( mHaveInterpolation.test( k ) )
                    {
                        delete mInterpolations( k );
                    }
                }

                // delete container
                mInterpolations.clear();
            }
        }

// ----------------------------------------------------------------------------

        bool has_interpolation(uint aBSplineMeshIndex)
        {
           return mHaveInterpolation.test(aBSplineMeshIndex);
        }

// ----------------------------------------------------------------------------

        void init_interpolation( uint aBSplineMeshIndex )
        {
            if ( ! mHaveInterpolationContainer )
            {
                mInterpolations.resize( gNumberOfMeshes, nullptr );    //FIXME make this size flexible to what is needed
                mHaveInterpolationContainer = true;
            }

            MORIS_ASSERT( aBSplineMeshIndex < gNumberOfMeshes, "Number of BSpline meshes on this Lagrange mesh exceeds gNumberOfMeshes. Increase gNumberOfMeshes" );
            MORIS_ASSERT( ! mHaveInterpolation.test( aBSplineMeshIndex ), "tried to intit an interpolation object that already exists" );    //FIXME Mathias

            mInterpolations( aBSplineMeshIndex ) = new Lagrange_Node_Interpolation;
            mHaveInterpolation.set( aBSplineMeshIndex );
        }

// ----------------------------------------------------------------------------

        /**
         * MTK Interface: return the coords of this node as Moris::Mat
         */
        Matrix< DDRMat > get_coords() const
        {
            Matrix< DDRMat > aCoords( 1, N );
            for( uint k=0; k<N; ++k )
            {
                aCoords( k ) = mXYZ[ k ];
            }
            return aCoords;
        }

// ----------------------------------------------------------------------------
        /**
         * Returns an array of size [N] telling the proc local ijk-position
         * of the node  on the current level.
         *
         * @return luint pointer to array containing ijk-position
         *               careful: node must not go out of scope.
         */
        const luint * get_ijk( ) const
        {
            return mIJK;
        }

// ----------------------------------------------------------------------------

        /**
         * set XYZ coordinates
         *
         * @param[in] aXYZ    array containing coordinates
         *
         * @return void
         */
         void set_xyz( const real * aXYZ )
         {
             // save ijk position in memory.
             for( uint k=0; k<N; ++k )
             {
                 mXYZ[ k ] = aXYZ[ k ];
             }
         }

// ----------------------------------------------------------------------------

         /**
          * Returns an array of size [N] telling the xyz-position
          * of the node
          *
          * @return double pointer to array containing xyz-position
          *               careful: node must not go out of scope.
          */
         const real * get_xyz() const
         {
             return mXYZ;
         }

// ----------------------------------------------------------------------------

         /**
          * set the T-Matrix coefficients
          */
         //void set_t_matrix( const Matrix< DDRMat > & aTMatrix )
         //{
         //    mTMatrix = aTMatrix;
         //}

// ----------------------------------------------------------------------------

         /**
          * set the DOFs
          */
         void set_coefficients( const uint                   aBSplineMeshIndex,
                 Vector< mtk::Vertex* > & aDOFs )
         {
             mInterpolations( aBSplineMeshIndex )->set_coefficients( aDOFs );
         }

// ----------------------------------------------------------------------------

         /**
          * set the weights
          */
         void set_weights( const uint               aBSplineMeshIndex,
                           const Matrix< DDRMat > & aWeights )
         {
             mInterpolations( aBSplineMeshIndex )->set_weights( aWeights );
         }

// ----------------------------------------------------------------------------

         /**
          * return a pointer to the interpolation object
          */
         mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
         {
             MORIS_ASSERT( mHaveInterpolation.test( aBSplineMeshIndex),
                      "tried to access an interpolation object for vertex ID %-5i and Index %-5i that does not exist"
                      , this->get_id(), this->get_index() );

             return mInterpolations( aBSplineMeshIndex );
         }

// ----------------------------------------------------------------------------

         /**
          * return a pointer to the interpolation object ( const version )
          */
         const mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const
         {
             MORIS_ASSERT( mHaveInterpolation.test( aBSplineMeshIndex),
                      "tried to access an interpolation object for vertex ID %-5i and Index %-5i that does not exist"
                      , this->get_id(), this->get_index() );

             return mInterpolations( aBSplineMeshIndex );
         }

// ----------------------------------------------------------------------------

         /**
          * return the T-Matrix coefficients
          */
         //const Matrix< DDRMat > * get_t_matrix() const
         //{
         //    return & mTMatrix;
         //}

// ----------------------------------------------------------------------------

         /**
          * return the DOF pointers
          */
         // const Cell< mtk::Vertex* > & get_adof_pointers() const
         //{
         //    return mDOFs;
         //}

// ----------------------------------------------------------------------------

         /**
          * return the DOF pointers
          */
         //Cell< mtk::Vertex* > & get_adof_pointers()
         //{
         //    return mDOFs;
         //}

// ----------------------------------------------------------------------------

         /**
          * return the IDs of used basis
          */
         //Mat< moris_id > get_adof_ids() const
         //{
         //    // allocate matrix with IDs
         //    uint tNumberOfDOFs = mDOFs.size();

         //    // create output matrix
         //    Mat< moris_id > aIDs( tNumberOfDOFs, 1 );
         //
         //    // write ids into matrix
         //    for( uint k=0; k<tNumberOfDOFs; ++k )
         //    {
         //        aIDs( k ) = mDOFs( k )->get_id();
         //    }
         //
         //    return aIDs;
        // }

// ----------------------------------------------------------------------------

         /**
          * return the indices of used basis
          */
         //Mat< moris_index > get_adof_indices() const
         //{
         //    // allocate matrix with IDs
         //    uint tNumberOfDOFs = mDOFs.size();

         //    // create output matrix
         //    Mat< moris_index > aIDs( tNumberOfDOFs, 1 );

         //    // write ids into matrix
         //    for( uint k=0; k<tNumberOfDOFs; ++k )
         //    {
         //        aIDs( k ) = mDOFs( k )->get_index();
         //    }

         //    return aIDs;
         //}

// ----------------------------------------------------------------------------

         /**
          * return the owners of used basis
          */
         // Matrix< DDUMat > get_adof_owners() const
         // {
         //    // allocate matrix with IDs
         //    uint tNumberOfDOFs = mDOFs.size();
         //
         // create output matrix
         // Matrix< DDUMat > aOwners( tNumberOfDOFs, 1 );
         //
         // write ids into matrix
         // for( uint k=0; k<tNumberOfDOFs; ++k )
         // {
         //     aOwners( k ) = mDOFs( k )->get_owner();
         // }
         //
         // return aOwners;
         //}

// ----------------------------------------------------------------------------
    };

// ----------------------------------------------------------------------------

} /* namespace moris */
#endif /* SRC_HMR_CL_HMR_LAGRANGE_NODE_HPP_ */

