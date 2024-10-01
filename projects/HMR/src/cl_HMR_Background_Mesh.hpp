/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Mesh.hpp
 *
 */

#pragma once

#include <string>

#include "cl_HMR_Background_Edge.hpp"            //HMR/src
#include "cl_HMR_Background_Element.hpp"         //HMR/src
#include "cl_HMR_Background_Element_Base.hpp"    //HMR/src
#include "cl_HMR_Background_Facet.hpp"           //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp"       //HMR/src
#include "cl_HMR_Domain.hpp"                     //HMR/src
#include "cl_HMR_Parameters.hpp"                 //HMR/src
#include "HMR_Globals.hpp"                       //HMR/src
#include "HMR_Tools.hpp"                         //HMR/src
#include "assert.hpp"
#include "cl_Communication_Tools.hpp"            //COM/src
#include "cl_Communication_Manager.hpp"          //COM/src
#include "cl_Tracer.hpp"
#include "fn_stringify_matrix.hpp"

#include "fn_equal_to.hpp"     //COM/src

#include "moris_typedefs.hpp"        //COR/src
#include "cl_Vector.hpp"             //CNT/src

#include "cl_Stopwatch.hpp"    //CHR/src

#include "cl_Matrix.hpp"       //LINALG/src

namespace moris::hmr
{
    //--------------------------------------------------------------------------------

    /**
     * \brief Background Mesh Class which is templated against dimension.
     *        To be created by the factory.
     */
    template< uint N >
    class Background_Mesh : public Background_Mesh_Base
    {
        //! contains ijk lookup tables for the whole mesh
        const Domain< N > mDomain;

        //! contains ijk lookup tables for the proc local subdomain
        Domain< N > mMySubDomain;

        //! contains the length of an element
        real mElementLength[ gMaxNumberOfLevels ][ N ];

        //! contains the coordinates of first node on aura
        real mDomainOffset[ N ];

        //--------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------
        /**
         * constructor for templated mesh class
         *
         * param[in]  aParameters   pointer to user defined settings container
         */
        Background_Mesh( const Parameters* aParameters )
                : Background_Mesh_Base( aParameters )
                , mDomain( aParameters->get_domain_ijk(), mPaddingSize )
        {
            // Log Trace this function
            Tracer tTracer( "HMR", "Background Mesh", "Create" );

            // create mesh decomposition
            this->decompose_mesh();

            // calculate element with, length and height for each level
            this->calculate_element_length();

            // calculate coordinate of first node on proc ( is in aura )
            this->calculate_domain_offset();

            // initialize first layer of elements
            this->initialize_coarsest_elements();

            // indices for elements on coarsest level within proc domain
            this->create_coarsest_frame();

            // set element properties on coarsest level // create aura and inverse aura
            // assumption: aura size = padding size
            this->finalize_coarsest_elements();

            // synchronize with other procs
            this->synchronize_coarsest_aura();

            // populates padding element container
            this->collect_coarsest_padding_elements();

            // create list of active elements
            this->collect_active_elements();

            // create list of active elements including aura
            this->collect_active_elements_including_aura();

            // calculate neighborhood
            this->collect_neighbors();

            // calculate indices for elements
            this->update_element_indices();
        }

        //--------------------------------------------------------------------------------

        /**
         * Mesh destructor. Needs to destroy all element pointers on coarsest
         * level. Higher level elements are destroyed implicitly by element
         * destructor.
         */
        ~Background_Mesh() override
        {
            // delete pointers in element cell
            for ( auto p : mCoarsestElementsIncludingAura )
            {
                delete p;
            }
        }

        //--------------------------------------------------------------------------------

        /**
         * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
         *                                       * < max number of levels >
         *
         * @return         Matrix< DDLUMat > number of elements per direction on
         *                              proc, including aura
         */
        Matrix< DDLUMat >
        get_number_of_elements_per_direction_on_proc() const override
        {
            Matrix< DDLUMat > aMat( N, gMaxNumberOfLevels );

            for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < N; ++k )
                {
                    aMat( k, l ) = mMySubDomain.mNumberOfElementsPerDimension[ l ][ k ];
                }
            }
            return aMat;
        }

        //--------------------------------------------------------------------------------

        /**
         * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
         *                                            * < max number of levels >
         *
         * @return         Matrix< DDLUMat > number of elements per direction
         *                                   within whole mesh, including aura
         */
        Matrix< DDLUMat >
        get_number_of_elements_per_direction() const override
        {
            Matrix< DDLUMat > aMat( N, gMaxNumberOfLevels );

            for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < N; ++k )
                {
                    aMat( k, l ) = mDomain.mNumberOfElementsPerDimension[ l ][ k ];
                }
            }
            return aMat;
        }

        //--------------------------------------------------------------------------------

        /**
         * Returns a Matrix< DDLUMat > containing the ijk positions of the calculation
         *                        domain on the proc
         *
         * @return Matrix< DDLUMat >
         */
        Matrix< DDLUMat >
        get_subdomain_ijk() const override
        {
            Matrix< DDLUMat > aMat( 2, N );

            for ( uint k = 0; k < N; ++k )
            {
                aMat( 0, k ) = mMySubDomain.mFrameIJK[ 0 ][ k ][ 0 ];
                aMat( 1, k ) = mMySubDomain.mFrameIJK[ 0 ][ k ][ 1 ];
            }

            return aMat;
        }

        //--------------------------------------------------------------------------------

        /**
         * Returns the ijk-offset of domain of current proc.
         * This value is needed to transform global IDs to local ones and
         * vice versa
         *
         * @return Matrix< DDLUMat > of dimension  < number of dimensions >
         *                                  * <max number of levels>
         */
        Matrix< DDLUMat >
        get_subdomain_offset_of_proc() override
        {
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            Matrix< DDLUMat > aIJK( tNumberOfDimensions, gMaxNumberOfLevels );

            for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < tNumberOfDimensions; ++k )
                {
                    aIJK( k, l ) = mMySubDomain.mAuraIJK[ l ][ k ][ 0 ];
                }
            }

            return aIJK;
        }
        //--------------------------------------------------------------------------------

        /**
         * calculates the node coordinates of an element
         *
         * @param[in]   aElement    Element to be processed
         * @param[out]  aNodeCoords Matrix containing the node coordinates
         *                          ( 3 x 2^n )
         */
        void calc_corner_nodes_of_element(
                const Background_Element_Base* aElement,
                Matrix< DDRMat >&              aNodeCoords ) override;

        //--------------------------------------------------------------------------------

        /**
         * calculates the coordinates of the center of the element
         *
         * @param[in]   aElement    Element to be processed
         * @param[out]  aNodeCoords Matrix containing the node coordinates
         *                          ( 3 x 1 )
         */
        void calc_center_of_element(
                const Background_Element_Base* aElement,
                Matrix< DDRMat >&              aNodeCoords ) override;

        //--------------------------------------------------------------------------------

        /**
         * returns a pointer to an element on level Zero
         *
         * @param[in] aI proc local i-position of element
         * @param[in] aJ proc local j-position of element
         *
         * @return pointer to element on top level
         *
         */
        Background_Element_Base*
        get_coarsest_element_by_ij(
                luint aI,
                luint aJ ) override
        {
            return mCoarsestElementsIncludingAura( this->calc_subdomain_id_of_element( 0, aI, aJ ) );
        }

        //--------------------------------------------------------------------------------

        /**
         * returns a pointer to an element on level Zero
         *
         * @param[in] aI proc local i-position of element
         * @param[in] aJ proc local j-position of element
         *
         * @return pointer to element on top level
         *
         */
        Background_Element_Base*
        get_coarsest_element_by_ijk(
                luint aI,
                luint aJ,
                luint aK ) override
        {
            return mCoarsestElementsIncludingAura( this->calc_subdomain_id_of_element( 0, aI, aJ, aK ) );
        }

        //--------------------------------------------------------------------------------

        /**
         * returns the offset of the current proc
         *
         * @return Matrix< DDRMat >
         */
        Matrix< DDRMat >
        get_domain_offset() override
        {
            Matrix< DDRMat > aMat( N, 1 );

            for ( uint k = 0; k < N; ++k )
            {
                aMat( k ) = mDomainOffset[ k ];
            }

            return aMat;
        }

        //--------------------------------------------------------------------------------

      protected:
        //--------------------------------------------------------------------------------

        /**
         * Calculates a local ID from a given level
         * and global ID
         *
         * @param[in]    aLevel   level of element to be investigated
         * @param[in]    aID      global ID of element to be investigated
         *
         * @return       luint    local ID on submesh of proc
         */
        luint calc_subdomain_id_from_global_id(
                uint  aLevel,
                luint aID ) const override;

        //--------------------------------------------------------------------------------

        /**
         * @brief This function takes mDomainID as input and computes the global ijk of the element
         *
         * @param aLevel
         * @param aID
         * @param aIJK
         */

        void calc_ijk_from_global_id(
                const uint&  aLevel,
                const luint& aID,
                luint*       aIJK ) const override;

        //--------------------------------------------------------------------------------

        /**
         * Unets element active flag, sets element flag,
         * initializes children, and removes refinement queue flag
         *
         * param[in]    aElement  element to be refined
         * param[in]    aKeepState  Indicates
         *
         * @return      vpid
         */
        void refine_element( Background_Element_Base* aElement, const bool aKeepState ) override;

        //--------------------------------------------------------------------------------

        /**
         *
         * Protected function that checks neighbors and inserts them into
         * aElement. This function is to be used for elements on level zero
         * only.
         *
         * @param[inout]  aElement    element to be processed
         *
         * @return void
         *
         */
        void collect_neighbors_on_level_zero() override;

        //--------------------------------------------------------------------------------

        /**
         * calculate global ID from level and from i-position ( 1D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         *
         * @return luint global ID of element
         *
         */
        luint calc_domain_id_of_element(
                uint  aLevel,
                luint aI ) const override;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculate global ID from level and from ij-position ( 2D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         * @param[in] aJ     j-index of element
         *
         * @return luint global ID of element
         *
         */
        luint calc_domain_id_of_element(
                uint  aLevel,
                luint aI,
                luint aJ ) const override;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculate global ID from level and from ijk-position ( 3D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         * @param[in] aJ     j-index of element
         * @param[in] aK     k-index of element
         *
         * @return luint global ID of element
         *
         */
        luint calc_domain_id_of_element(
                uint  aLevel,
                luint aI,
                luint aJ,
                luint aK ) const override;
        //--------------------------------------------------------------------------------

        /**
         * subroutine for collect_side_set that collects elements on coarsest
         * level for a side
         */
        void collect_coarsest_elements_on_side(
                uint                              aSideOrdinal,
                Vector< Background_Element_Base* >& aCoarsestElementsOnSide ) override;

        //--------------------------------------------------------------------------------

      private:
        //--------------------------------------------------------------------------------

        /**
         * calculates the element length lookup table
         *
         * @return void
         */
        void
        calculate_element_length()
        {
            // get domain dimensions from settings
            const Vector< real >& tDomainDimensions = mParameters->get_domain_dimensions();

            // get number of elements on coarsest level from settings
            const Vector< uint >& tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

            // calculate width for first level
            Matrix< DDRMat > tInitElemSize( 3, 1, 0.0 );
            for ( uint iDim = 0; iDim < N; ++iDim )
            {
                mElementLength[ 0 ][ iDim ] = tDomainDimensions( iDim ) / ( (real)( tNumberOfElements( iDim ) ) );
                tInitElemSize( iDim )       = mElementLength[ 0 ][ iDim ];
            }

            // print the coarsest element size to screen
            MORIS_LOG_SPEC( "Initial BG element size", ios::stringify_log( tInitElemSize ) );

            // loop over all higher levels
            for ( uint iLevel = 1; iLevel < gMaxNumberOfLevels; ++iLevel )
            {
                // loop over all dimensions
                for ( uint iDim = 0; iDim < N; ++iDim )
                {
                    // calculate length of element
                    mElementLength[ iLevel ][ iDim ] = 0.5 * mElementLength[ iLevel - 1 ][ iDim ];
                }
            }
        }

        //--------------------------------------------------------------------------------

        /**
         * calculates the coordinate of first node on proc
         * ( including aura )
         *
         * @return void
         */
        void
        calculate_domain_offset()
        {
            // get domain offset
            Vector< real > tParametersOffset = mParameters->get_domain_offset();

            // get padding size
            real tPaddingSize = (real)mParameters->get_padding_size();

            // subtract padding size from offset
            for ( uint k = 0; k < N; ++k )
            {
                mDomainOffset[ k ] = tParametersOffset( k ) - tPaddingSize * mElementLength[ 0 ][ k ];
            }
        }

        //--------------------------------------------------------------------------------

        /**
         * Processor decomposition method to minimize the mesh interface between adjacent
         * processors. This functions is called for minimum mesh interface decomposition
         * method only (mProcDecompMethod = 2 in hmr::Parameters).
         *
         * @param[in]  aNumberOfDimensions   must be 1, 2, or 3
         * @param[in]  aMeshDims             the dimensions of the mesh
         * @param[out] aProcDims             determined processor grid dimensions
         *
         * @return void
         */
        void
        create_proc_dims_min_mesh_interface(
                uint            aNumberOfDimensions,
                Vector< uint >& aMeshDims,
                Vector< uint >& aProcDimensions )
        {

            // This function determines the processor dimensions based on minimizing mesh grid interfaces.
            // Eg. If the mesh is 100x10 elements, and we have 10 processors, the processor layout will be 10x1.

            real tInterfaceCount = 0;

            // 1D Processor Grid
            if ( aNumberOfDimensions == 1 )
            {
                aProcDimensions.resize( 1 );
                aProcDimensions( 0 ) = par_size();
            }

            // 2D Processor Grid
            else if ( aNumberOfDimensions == 2 )
            {
                aProcDimensions.resize( 2 );
                // Iterating through "i" (x) direction processors
                for ( uint i = 1; i <= (uint)par_size(); i++ )
                {

                    // initializing interface count value
                    if ( i == 1 )
                    {
                        tInterfaceCount = aMeshDims( 1 ) * ( i - 1 ) + aMeshDims( 0 ) * ( par_size() / i - 1 );
                        // Assigning processor grid dimensions
                        aProcDimensions( 0 ) = i;
                        aProcDimensions( 1 ) = par_size() / i;
                    }

                    // Otherwise, is i a factor of the number of processors?
                    else if ( (uint)par_size() % i == 0 )
                    {
                        // Is the last calculated interface count larger than the current iterations?
                        if ( tInterfaceCount > ( aMeshDims( 1 ) * ( i - 1 ) + aMeshDims( 0 ) * ( par_size() / i - 1 ) ) )
                        {
                            tInterfaceCount = aMeshDims( 1 ) * ( i - 1 ) + aMeshDims( 0 ) * ( par_size() / i - 1 );
                            // Assigning processor grid dimensions
                            aProcDimensions( 0 ) = i;
                            aProcDimensions( 1 ) = par_size() / i;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }

            // 3D Processor Grid
            else if ( aNumberOfDimensions == 3 )
            {
                aProcDimensions.resize( 3 );

                // Iterating through i-direction processor possibilities
                for ( uint i = 1; i <= (uint)par_size(); i++ )
                {
                    // Is "i" a factor of total processors used?
                    if ( (uint)par_size() % i == 0 )
                    {
                        // Iterating through j-direction processors possibilities
                        for ( uint j = 1; j <= (uint)par_size() / i; j++ )
                        {
                            // Initializing interface count value
                            if ( i == 1 && j == 1 )
                            {
                                tInterfaceCount = aMeshDims( 0 ) * aMeshDims( 1 ) * ( par_size() / ( i * j ) - 1 ) +    //
                                                  aMeshDims( 0 ) * aMeshDims( 2 ) * ( j - 1 ) + aMeshDims( 1 ) * aMeshDims( 2 ) * ( i - 1 );

                                // Replace processor grid dimensions
                                aProcDimensions( 0 ) = i;
                                aProcDimensions( 1 ) = j;
                                aProcDimensions( 2 ) = par_size() / ( i * j );
                            }

                            // Is i*j a factor of processor count?
                            else if ( (uint)par_size() % ( i * j ) == 0 )
                            {

                                // If the old interface count is higher than the current iteration
                                if ( tInterfaceCount >                                                          //
                                        ( aMeshDims( 0 ) * aMeshDims( 1 ) * ( par_size() / ( i * j ) - 1 ) +    //
                                                aMeshDims( 0 ) * aMeshDims( 2 ) * ( j - 1 ) + aMeshDims( 1 ) * aMeshDims( 2 ) * ( i - 1 ) ) )
                                {

                                    // Replace interface count
                                    tInterfaceCount = aMeshDims( 0 ) * aMeshDims( 1 ) * ( par_size() / ( i * j ) - 1 ) +    //
                                                      aMeshDims( 0 ) * aMeshDims( 2 ) * ( j - 1 ) + aMeshDims( 1 ) * aMeshDims( 2 ) * ( i - 1 );

                                    // Repace processor grid dimensions
                                    aProcDimensions( 0 ) = i;
                                    aProcDimensions( 1 ) = j;
                                    aProcDimensions( 2 ) = par_size() / ( i * j );
                                }
                            }
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------
        /**
         * Calculates ijk lookup table for current proc.
         *
         * @return  void
         */
        void
        decompose_mesh()
        {
            // print output info
            MORIS_LOG_INFO( "Decomposing the mesh over %u processor(s).", (unsigned int)par_size() );

            // Processor decomposition method 0=UserDefined, 1=MPI (original) 2=min mesh interface 3=manual
            uint tDecompMethod = mParameters->get_processor_decomp_method();

            // Pulling mesh dimensions from parameters
            Vector< uint > tNumberOfElementsPerDimension = mParameters->get_number_of_elements_per_dimension();

            switch ( tDecompMethod )
            {
                // User defined processor grid
                case 0:
                {
                    mProcessorDimensions = mParameters->get_processor_dimensions();

                    // Checking if user defined processor dimensions matches mesh dimensions, N.
                    MORIS_ERROR(
                            mProcessorDimensions.max() == N,
                            "hmr::Background_Mesh::decompose_mesh() - "
                            "User defined processor grid dimensions incompatible with mesh dimensions." );

                    // Calculating the product of user defined proc dims dimensions
                    uint tProcCount = 1;
                    for ( uint iDim = 0; iDim < N; ++iDim )    // N is number of spatial dimensions
                    {
                        tProcCount = tProcCount * mProcessorDimensions( iDim );
                    }

                    // check that the total number of processors for user-defined mesh splitting matches with actual number of processors
                    if ( (uint)par_size() != tProcCount )
                    {
                        MORIS_ERROR(
                                false,
                                "hmr::Background_Mesh::decompose_mesh() - "
                                "User defined processor grid dimensions do not match number of processors used." );
                    }

                    break;
                }

                // MPI decomposition method. Minimize processor interface
                case 1:
                {
                    // No alteration needed.  Processor decomposition occurs within create_proc_cart
                    break;
                }

                // Minimize mesh interface decomposition method
                case 2:
                {
                    create_proc_dims_min_mesh_interface(
                            N,
                            tNumberOfElementsPerDimension,
                            mProcessorDimensions );
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "hmr::Background_Mesh::decompose_mesh() - Invalid processor decomposition method defined." );
                    break;
                }
            }    // end switch: mesh decomposition method

            create_proc_cart(
                    tDecompMethod,
                    N,
                    mProcessorDimensions,
                    mMyProcCoords,
                    mMyProcNeighbors );

            // reporting the proc dims to the logger
            if ( par_rank() == 0 )
            {
                // how many dimensions?
                switch ( N )
                {
                    case 2:
                    {
                        MORIS_LOG_INFO( "Background mesh split (x,y) across processors is (%i,%i)",
                                mProcessorDimensions( 0 ),
                                mProcessorDimensions( 1 ) );
                        break;
                    }
                    case 3:
                    {
                        MORIS_LOG_INFO( "Background mesh split (x,y,z) across processors is (%i,%i,%i)",
                                mProcessorDimensions( 0 ),
                                mProcessorDimensions( 1 ),
                                mProcessorDimensions( 2 ) );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "hmr::Background_Mesh::decompose_mesh() - Invalid number of spatial dimensions." );
                    }
                }
            }

            // calculate number of elements per dimension
            Vector< uint > tNumberOfElementsPerDimensionOnProc( N );

            // remainder used if tNumberOfElementsPerDimension(k) isn't a
            // multiple of mProcDims(k).
            Vector< uint > tRemainder( N, 1 );

            for ( uint k = 0; k < N; ++k )
            {
                /*
                 * This conditional determines whether a remainder element
                 * is added to the processor dimension.  This will only
                 * affect situations where the number of elements in a
                 * direction is not a multiple of the number of processors
                 * in that direction.
                 *
                 * For example, if there are 10 elements and 3 processors
                 * in a direction, this will assign 4 elements to proc 1
                 * and 3 elements to the other 2.
                 */
                if ( mMyProcCoords( k ) < ( tNumberOfElementsPerDimension( k ) % mProcessorDimensions( k ) ) )
                {
                    // Remainder applied
                    tRemainder( k ) = 1;
                }
                else
                {
                    // Remainder omitted
                    tRemainder( k ) = 0;
                }

                // assigning number of elements on this direction of the processor
                tNumberOfElementsPerDimensionOnProc( k ) =
                        ( tNumberOfElementsPerDimension( k ) - ( tNumberOfElementsPerDimension( k ) % mProcessorDimensions( k ) ) ) / mProcessorDimensions( k ) + tRemainder( k );
            }

            // calculate decomposition domain
            // set owned and shared limits
            Matrix< DDLUMat > tDomainIJK( 2, N );

            sint tNumElemsOnProc   = 1;
            sint tNumTotalElements = 1;

            for ( uint k = 0; k < N; ++k )
            {
                // calculates domain start taking into account remainder elements
                tDomainIJK( 0, k ) =
                        mPaddingSize + tNumberOfElementsPerDimensionOnProc( k ) * mMyProcCoords( k ) +    //
                        ( 1 - tRemainder( k ) ) * ( tNumberOfElementsPerDimension( k ) % mProcessorDimensions( k ) );

                tDomainIJK( 1, k ) = tDomainIJK( 0, k ) + tNumberOfElementsPerDimensionOnProc( k ) - 1;

                // compute total number of elements (based on input specs)
                tNumTotalElements *= tNumberOfElementsPerDimension( k );

                // compute number of elements on proc
                tNumElemsOnProc *= tNumberOfElementsPerDimensionOnProc( k );
            }

            // check that total number of elements across all processor equals required number
            sint tNumElemsOnAllProcs = sum_all( tNumElemsOnProc );

            MORIS_ERROR(
                    tNumElemsOnAllProcs == tNumTotalElements,
                    "Total number of elements on all processors: %d (should be %d).\n",
                    tNumElemsOnAllProcs,
                    tNumTotalElements );

            // create proc domain
            Domain< N > tMySubDomain( tDomainIJK, mPaddingSize );

            // move domain object to member
            mMySubDomain = std::move( tMySubDomain );

            // print proc area
            if ( gLogger.get_severity_level() < 1 )
            {
                moris_id tMyRank = par_rank();

                uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

                std::string tString = "Processor #" + std::to_string( tMyRank );

                // add dots to the string for pretty output
                if ( tMyRank < 10 )
                {
                    tString += " ... :";
                }
                else if ( tMyRank < 100 )
                {
                    tString += " .. :";
                }
                else if ( tMyRank < 1000 )
                {
                    tString += " . :";
                }
                else if ( tMyRank < 10000 )
                {
                    tString += "  :";
                }
                else
                {
                    tString += " :";
                }

                if ( tNumberOfDimensions == 1 )
                {
                    // MORIS_LOG_INFO("%s owns i domain ", tString.c_str() ) ;
                    tString += " owns i domain";
                }
                else if ( tNumberOfDimensions == 2 )
                {
                    // MORIS_LOG_INFO("%s owns i-j domain ", tString.c_str()  ) ;
                    tString += " owns i-j domain";
                }
                else if ( tNumberOfDimensions == 3 )
                {
                    // MORIS_LOG_INFO("%s owns i-j-k domain ", tString.c_str() ) ;
                    tString += " owns i-j-k domain";
                }

                // print ijk domain
                // MORIS_LOG_INFO("%lu-%lu",
                //         ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 0 ],
                //         ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 1 ] );

                std::string tLogString = std::to_string( (long unsigned int)mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 0 ] );
                tLogString += "-" + std::to_string( (long unsigned int)mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 1 ] );

                for ( uint iDim = 1; iDim < N; ++iDim )
                {
                    // MORIS_LOG_INFO(", %lu-%lu ",
                    //         ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ iDim ][ 0 ],
                    //         ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ iDim ][ 1 ] );
                    tLogString += ", " + std::to_string( (long unsigned int)mMySubDomain.mDomainIJK[ 0 ][ iDim ][ 0 ] );
                    tLogString += "-" + std::to_string( (long unsigned int)mMySubDomain.mDomainIJK[ 0 ][ iDim ][ 1 ] );
                }
                MORIS_LOG_INFO_ALL_PROCS( "%s: %s", tString.c_str(), tLogString.c_str() );
            }

            // test if settings are OK
            if ( par_size() != 1 )
            {
                // get number of dimensions from settings
                uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

                Matrix< DDLUMat > tProcSplit( tNumberOfDimensions, 1 );

                int tError = 0;

                // loop over all dimensions and copy elements
                for ( uint k = 0; k < tNumberOfDimensions; ++k )
                {
                    // do not consider pseudo 1D directions or 2 Procs in 1 direction since aura overlap isn't an issue
                    if ( mProcessorDimensions( k ) > 2 )
                    {
                        // determine how many elements are on this proc in this domain minus the padding elements in that direction
                        tProcSplit( k ) = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ k ] - 2 * mPaddingSize;

                        // there should be at least on element left
                        MORIS_ERROR(
                                tProcSplit( k ) > 0,
                                "Processor %u does not have any non-padding/aura elements in direction %u. \n",
                                par_rank(),
                                k );

                        if ( mParameters->use_number_aura() )
                        {
                            // the aura of next-adjacent procs cannot use the same node since they do not communicate
                            if ( tProcSplit( k ) < ( mPaddingSize * 2 ) + 1 )
                            {
                                tError = ( mPaddingSize * 2 ) + 1 - tProcSplit( k );
                            }
                        }

                        // if not using numbered aura, then the remainder just has to be larger than the padding/aura
                        else
                        {
                            if ( tProcSplit( k ) < mPaddingSize )
                            {
                                tError = mPaddingSize - tProcSplit( k );
                            }
                        }
                    }
                }

                // get error code across all processors
                int tTotalError = max_all( tError );

                // test if error occurred
                if ( tTotalError > 0 )
                {
                    for ( uint k = 0; k < tNumberOfDimensions; ++k )
                    {
                        if ( tProcSplit( k ) < mPaddingSize )
                        {
                            // print in correct grammar
                            if ( tProcSplit( k ) == 1 )
                            {
                                MORIS_LOG_INFO( "In dimension %u, each proc domain gets 1 element, need at least %lu.\n",
                                        (unsigned int)k,
                                        (long unsigned int)mPaddingSize );
                            }
                            else
                            {
                                MORIS_LOG_INFO( "In dimension %u, each proc domain gets %lu elements, need at least %lu.\n",
                                        (unsigned int)k,
                                        (long unsigned int)tProcSplit( k ),
                                        (long unsigned int)mPaddingSize );
                            }
                        }
                    }

                    barrier();

                    MORIS_ERROR(
                            false,
                            "Mesh too coarse for selected split and padding size; check individual processors; maximum mismatch: %d.\n",
                            tTotalError );
                }
            }
        }    // end function: hmr::Background_Mesh::decompose_mesh()

        //--------------------------------------------------------------------------------
        /**
         * In this subroutine, the coarsest layer of elements on the proc,
         * including the aura, is generated. The element pointers are created
         * and stored in a fixed size cell, mCoarsestElementsIncludingAura.
         * New elements are created with the assumption that they are active,
         * the ID of the proc owning that element is not known yet, and set to
         * MORIS_UINT_MAX.
         *
         * @return  void
         */
        void initialize_coarsest_elements();

        //--------------------------------------------------------------------------------
        /**
         * This function loops over all elements on the coarsest level,
         * identifies padding elements and calculates element ownership. It
         * also identifies the elements on the coarsest level that belong to
         * the aura. In fact, there are two auras. The normal aura around the
         * proc domain, which contains the elements that belong to the neighbor
         * proc and are shared with the current proc, and the inverse aura,
         * which contains the elements that belong to the current proc,
         * and are shared with the neighbor. In most cases, padding elements are
         * not considered for the aura.
         *
         * @return  void
         */
        void finalize_coarsest_elements();

        //--------------------------------------------------------------------------------
        /**
         * Sometimes, it is more helpful to only access the elements that are
         * within the calculation domain of the proc, excluding the aura.
         * This function creates a list of local indices of these elements.
         *
         * @return  void
         */
        void create_coarsest_frame();

        //--------------------------------------------------------------------------------

        /**
         * Calculates number of elements on first level.
         * This function is only needed once.
         * Therefore, template specialization is not required.
         *
         * @return   Matrix< DDLUMat > of dimension < number of dimensions>
         *                       containing number of elements per direction
         *                       on coarsest proc, including aura
         */
        Matrix< DDLUMat >
        get_number_of_subdomain_elements_per_direction_on_level_zero()
        {
            Matrix< DDLUMat > aNumberOfElements( N, 1 );
            for ( uint k = 0; k < N; ++k )
            {
                aNumberOfElements( k ) = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ k ];
            }
            return aNumberOfElements;
        }

        //--------------------------------------------------------------------------------

        /**
         * Internal function called from initialize_coarsest_elements()/
         * Takes the pointer of the element and puts it into the correct
         * position in mCoarsestElementsIncludingAura
         *
         *
         */
        void
        insert_zero_level_element( 
                luint aPosition,
                Background_Element_Base* aElement )
        {
            mCoarsestElementsIncludingAura( aPosition ) = aElement;
        }

        //-------------------------------------------------------------------------------

        /**
         * calculate subdomain ID from level and from i-position ( 1D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         *
         * @return luint subdomain ID of element
         *
         */
        luint calc_subdomain_id_of_element(
                uint  aLevel,
                luint aI ) const;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculate subdomain ID from level and from ij-position ( 2D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         * @param[in] aJ     j-index of element
         *
         * @return luint subdomain ID of element
         *
         */
        luint calc_subdomain_id_of_element(
                uint  aLevel,
                luint aI,
                luint aJ ) const;
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculate subdomain ID from level and from ijk-position ( 3D case )
         *
         * @param[in] aLevel level of element
         * @param[in] aI     i-index of element
         * @param[in] aJ     j-index of element
         * @param[in] aK     k-index of element
         *
         * @return luint subdomain ID of element
         *
         */
        luint calc_subdomain_id_of_element(
                uint  aLevel,
                luint aI,
                luint aJ,
                luint aK ) const;

        //--------------------------------------------------------------------------------
        /**
         * Internal function called during element refinement.
         * Calculates IDs and Subdomain IDs for given ijk-positions
         *
         * @param[in]  aLevel         level of child elements to be considered
         * @param[in]  aIJK           matrix of dimension <number of dimensions>
         *                           * <number of children>
         * @param[out] aIDs           global IDs for each child element
         * @param[out] aSubdomainIDs  proc local IDs for each child element
         *
         * @return void
         */
        void calc_element_ids( uint      aLevel,
                const Matrix< DDLUMat >& aIJK,
                Matrix< DDLUMat >&       aIDs ) const;

        //--------------------------------------------------------------------------------

        void
        check_queued_element_for_padding( Background_Element_Base* aElement ) override
        {
            // only do something if this element belongs to me
            if ( aElement->get_owner() == mMyRank )
            {
                // get local ijk position of element
                const luint* tIJK = aElement->get_ijk();

                // get level of element
                uint tLevel = aElement->get_level();

                // perform aura check
                bool tIsPaddingCandidate = false;

                // loop over all dimensions
                for ( uint iDim = 0; iDim < N; ++iDim )
                {
                    // compute global coordinate of element
                    luint tI = tIJK[ iDim ] + mMySubDomain.mAuraIJK[ tLevel ][ iDim ][ 0 ];

                    // test if element is candidate for padding test
                    tIsPaddingCandidate = tIsPaddingCandidate
                                       || ( ( tI < mDomain.mDomainIJK[ tLevel ][ iDim ][ 0 ] + mPaddingRefinement )
                                               || ( tI > mDomain.mDomainIJK[ tLevel ][ iDim ][ 1 ] - mPaddingRefinement ) );
                }

                if ( tIsPaddingCandidate )
                {
                    // get neighbors from same level
                    Vector< Background_Element_Base* > tNeighbors;
                    aElement->get_neighbors_from_same_level( mPaddingRefinement, tNeighbors );

                    // loop over all neighbors
                    for ( auto tNeighbor : tNeighbors )
                    {
                        // test if neighbor is padding
                        if ( tNeighbor->is_padding() && !tNeighbor->is_queued_for_refinement() )
                        {
                            // flag padding element for refinement
                            tNeighbor->put_on_refinement_queue();
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------

        void
        get_element_in_bounding_box_memory_index(
                uint                           aPattern,
                const moris::Matrix< DDRMat >& aPoint,
                const moris::Matrix< DDRMat >& aBoundingBoxSize,
                moris::Matrix< DDLUMat >&      aElementMemoryIndex ) override
        {
            MORIS_ASSERT( par_size() == 1, "get_element_in_bounding_box_memory_index(), not tested in parallel" );
            uint tLevel = 0;
            //                for( uint Ik = 0; Ik < gMaxNumberOfLevels; Ik++ )
            //                {
            //                    if( equal_to( aBackgroundElementedgeLength( 0, 0 ), mElementLength[ Ik ][ 0 ] ) )
            //                    {
            //                        break;
            //                    }
            //
            //                    tLevel += 1;
            //                }

            // get domain offset
            Vector< real > tParametersOffset = mParameters->get_domain_offset();

            bool tCheck = true;

            sint tIJK[ N ];

            luint tBoundingBoxSizeIJK[ N ];
            luint tBoundingBoxStartEndIJK[ N ][ 2 ];
            for ( uint Ik = 0; Ik < N; ++Ik )
            {
                tIJK[ Ik ] = std::floor( ( aPoint( Ik ) - mDomainOffset[ Ik ] ) / mElementLength[ tLevel ][ Ik ] );

                if ( tIJK[ Ik ] < (sint)mDomain.mDomainIJK[ tLevel ][ Ik ][ 0 ] || tIJK[ Ik ] > (sint)mDomain.mDomainIJK[ tLevel ][ Ik ][ 1 ] )
                {
                    tCheck = false;
                    break;
                }

                tBoundingBoxSizeIJK[ Ik ] = std::ceil( aBoundingBoxSize( Ik ) / mElementLength[ tLevel ][ Ik ] );

                sint tStart = tIJK[ Ik ] - tBoundingBoxSizeIJK[ Ik ];

                if ( tStart >= (sint)mDomain.mDomainIJK[ tLevel ][ Ik ][ 0 ] )
                {
                    tBoundingBoxStartEndIJK[ Ik ][ 0 ] = tStart;
                }
                else
                {
                    tBoundingBoxStartEndIJK[ Ik ][ 0 ] = mDomain.mDomainIJK[ tLevel ][ Ik ][ 0 ];
                }

                sint tEnd = tIJK[ Ik ] + tBoundingBoxSizeIJK[ Ik ];

                if ( tEnd <= (sint)mDomain.mDomainIJK[ tLevel ][ Ik ][ 1 ] )
                {
                    tBoundingBoxStartEndIJK[ Ik ][ 1 ] = tEnd;
                }
                else
                {
                    tBoundingBoxStartEndIJK[ Ik ][ 1 ] = mDomain.mDomainIJK[ tLevel ][ Ik ][ 1 ];
                }
            }

            if ( tCheck )
            {
                Vector< Background_Element_Base* > tBackgroundElements;

                this->collect_coarsest_elements_in_bounding_box( tBackgroundElements, tBoundingBoxStartEndIJK, tLevel );

                luint tCounter = 0;

                // loop over frame and count active descendants
                for ( luint Ik = 0; Ik < tBackgroundElements.size(); ++Ik )
                {
                    tBackgroundElements( Ik )->get_number_of_active_descendants( aPattern, tCounter );
                }

                Vector< Background_Element_Base* > tActiveElements( tCounter, nullptr );

                tCounter = 0;

                for ( luint Ik = 0; Ik < tBackgroundElements.size(); ++Ik )
                {
                    tBackgroundElements( Ik )->collect_active_descendants( aPattern, tActiveElements, tCounter );
                }

                aElementMemoryIndex.set_size( tActiveElements.size(), 1 );

                for ( luint Ik = 0; Ik < tActiveElements.size(); ++Ik )
                {
                    aElementMemoryIndex( Ik ) = tActiveElements( Ik )->get_memory_index();
                }
            }
            else
            {
                aElementMemoryIndex.set_size( 0, 0 );
            }
        };

        //--------------------------------------------------------------------------------

        void collect_coarsest_elements_in_bounding_box( Vector< Background_Element_Base* >& aBackgroundElements,
                luint                                                                     aBoundingBoxStartEndIJK[][ 2 ],
                uint                                                                      alevel );

        //--------------------------------------------------------------------------------
    }; /* Background_Mesh */

    //--------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::initialize_coarsest_elements()
    {
        MORIS_ERROR( false, "Do not know how initialize elements\n" );
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::finalize_coarsest_elements()
    {
        MORIS_ERROR( false, "Don't know how to finalize coarsest level." );
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI ) const
    {
        MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called." );
        return 0;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ ) const
    {
        MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called." );
        return 0;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template< uint N >
    inline luint
    Background_Mesh< N >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ,
            luint aK ) const
    {
        MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called." );
        return 0;
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_subdomain_id_of_element(
            uint  aLevel,
            luint aI ) const
    {
        MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called." );
        return 0;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_subdomain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ ) const
    {
        MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called." );
        return 0;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_subdomain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ,
            luint aK ) const
    {
        MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called." );
        return 0;
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::calc_element_ids(
            uint                     aLevel,
            const Matrix< DDLUMat >& aIJK,
            Matrix< DDLUMat >&       aIDs ) const
    {
        MORIS_ERROR( false, "Don't know how to calculate ids yet." );
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline luint
    Background_Mesh< N >::calc_subdomain_id_from_global_id(
            uint  aLevel,
            luint aID ) const
    {
        MORIS_ERROR( false, "Don't know how to calculate IDs yet." );
        return 0;
    }

    //--------------------------------------------------------------------------------
    template< uint N >
    inline void
    Background_Mesh< N >::calc_ijk_from_global_id(
            const uint&  aLevel,
            const luint& aID,
            luint*       aIJK ) const
    {
        MORIS_ERROR( false, "Don't know how to calculate IDs yet." );
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::refine_element( Background_Element_Base* aElement, const bool aKeepState )
    {
        MORIS_ERROR( false, "Don't know how to refine element." );
    }

    //--------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::create_coarsest_frame()
    {
        MORIS_ERROR( false, "Don't know how to create coarsest frame." );
    }

    //-------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::collect_neighbors_on_level_zero()
    {
        MORIS_ERROR( false, "Don't know how to collect_neighbors_on_level_zero." );
    }

    //-------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::calc_corner_nodes_of_element(
            const Background_Element_Base* aElement,
            Matrix< DDRMat >&              aNodeCoords )
    {
        MORIS_ERROR( false, "Do not know how calculate corner nodes\n" );
    }

    //-------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::calc_center_of_element(
            const Background_Element_Base* aElement,
            Matrix< DDRMat >&              aNodeCoords )
    {
        MORIS_ERROR( false, "Do not know how calculate center of element\n" );
    }

    //-------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::collect_coarsest_elements_on_side(
            uint                              aSideOrdinal,
            Vector< Background_Element_Base* >& aCoarsestElementsOnSide )
    {
        MORIS_ERROR( false, "Do not know how to collect coarsest elements on side \n" );
    }

    //-------------------------------------------------------------------------------

    template< uint N >
    inline void
    Background_Mesh< N >::collect_coarsest_elements_in_bounding_box(
            Vector< Background_Element_Base* >& aBackgroundElements,
            luint                                    aBoundingBoxStartEndIJK[][ 2 ],
            uint                                     alevel )
    {
        MORIS_ERROR( false, "Do not know how initialize elements\n" );
    }

    //--------------------------------------------------------------------------------
}    // namespace moris::hmr

#include "cl_HMR_Background_Mesh_2D.hpp"
#include "cl_HMR_Background_Mesh_3D.hpp"
