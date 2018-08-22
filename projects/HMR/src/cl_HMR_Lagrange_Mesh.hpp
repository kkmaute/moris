/*
 * cl_HMR_Lagrange_Mesh.hpp
 *
 *  Created on: May 15, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_

#include "typedefs.hpp" //COR/src
#include "cl_Stopwatch.hpp" //CHR/src
//#include "cl_Map.hpp" //CON/src
#include "HMR_Globals.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Background_Element_Base.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Element.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

    /**
     * \brief the Lagrange_Mesh class calculates Lagrange nodes for a given
     *  background mesh.
     *
     */
    template< uint N, uint P >
    class Lagrange_Mesh : public Lagrange_Mesh_Base
    {
        //! Lookup table containing offset for node IDs
        luint mNodeLevelOffset[ gMaxNumberOfLevels ];

        //! Lookup table containing number of elements per dimension for each level
        luint mNumberOfElementsPerDimensionIncludingAura[ gMaxNumberOfLevels ][ N ];

        //! Lookup table for node IDs
        luint mMySubdomainOffset[ gMaxNumberOfLevels ][ N ];

// ----------------------------------------------------------------------------
    public:
// ----------------------------------------------------------------------------

        /**
         * Constructor for Lagrange Mesh
         *
         * @param[in] aParameters       ref to container of user defined settings
         * @param[in] aBackgroundMesh pointer to background mesh
         *
         */
        Lagrange_Mesh( const Parameters       * aParameters,
                       Background_Mesh_Base * aBackgroundMesh,
                       BSpline_Mesh_Base    * aBSplineMesh ) :
                       Lagrange_Mesh_Base( aParameters, aBackgroundMesh, aBSplineMesh, P )
        {

            // ask background mesh for number of elements per ijk-direction
            this->get_number_of_elements_per_dimension();

            // calculate lookup table mNodeLevelOffset
            this->calculate_level_offset();

            // find out coordinate of first point on proc subdomain
            this->calculate_subdomain_offset();

            // calculate any value that can change after refinement
            this->update_mesh();
        }

// ----------------------------------------------------------------------------

        /**
         * Default destructor.
         */
        ~Lagrange_Mesh()
        {
           this->delete_pointers();
        }


// ----------------------------------------------------------------------------

        /**
         * Creates a Lagrange element and links it to corresponding element
         * on background mesh.
         *
         * @param[in] aElement  pointer to element on background mesh
         *
         * @return Element*  new Lagrange element
         */
        Element *
        create_element( Background_Element_Base* aElement );


// ----------------------------------------------------------------------------
    private:
// ----------------------------------------------------------------------------


        /**
         * calculates domain wide unique node ID (1D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                const uint  & aLevel,
                const luint & aI )
        {
            if( aLevel < gMaxNumberOfLevels && N == 1 )
            {
                return  aI
                        + mMySubdomainOffset[ aLevel ][ 0 ]
                        + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculates domain wide unique node ID (2D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                const uint  & aLevel,
                const luint & aI,
                const luint & aJ )
        {
            if( aLevel < gMaxNumberOfLevels && N == 2 )
            {
                return  aI + mMySubdomainOffset[ aLevel ][ 0 ]
                        + ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] )
                        *( P*mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 0 ] + 1)
                        + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculates domain wide unique node ID (1D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @param[in]  aK        proc local k-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                const uint  & aLevel,
                const luint & aI,
                const luint & aJ,
                const luint & aK )
        {
            if( aLevel < gMaxNumberOfLevels && N == 3 )
            {
                return  aI + mMySubdomainOffset[ aLevel ][ 0 ]
                        + ( P*mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 0 ] + 1)
                        * ( ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] )
                        + (aK + mMySubdomainOffset[ aLevel ][ 2 ] ) * ( P*
                          mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 1 ] + 1) )
                        + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

// ----------------------------------------------------------------------------

        /**
         * Internal function. Asks the background mesh for number of elements
         * per direction, stores result in  mAuraNumberOfElementsPerDimension
         *
         * @return void
         *
         */
        void
        get_number_of_elements_per_dimension()
        {
            // get elements per level from background mesh
            Mat< luint > tMat =
                mBackgroundMesh->get_number_of_elements_per_direction();

            // convert matrix to fixed size array
            for( uint l=0; l<gMaxNumberOfLevels; ++l )
            {
                for( uint k=0; k<N; ++k )
                {
                    mNumberOfElementsPerDimensionIncludingAura[ l ][ k ]
                        = tMat( k, l );
                }
            }

        }

// ----------------------------------------------------------------------------

        /**
         *  Private function, creates the mNodeLevelOffset lookup table.
         *
         *  @return void
         */
        void
        calculate_level_offset()
        {
            // calculate node level offset
            mNodeLevelOffset[ 0 ] = 0;

            for( uint l=1; l<gMaxNumberOfLevels; ++l )
            {
                // calculate number of nodes on this level
                luint tNumberOfNodes = 1;
                for( uint k=0; k<N; ++k )
                {
                    tNumberOfNodes *= P*
                        mNumberOfElementsPerDimensionIncludingAura[ l-1 ][ k ] + 1;
                }

                // add number of nodes to offset table
                mNodeLevelOffset[ l ] =
                        mNodeLevelOffset[ l-1 ] + tNumberOfNodes;
            }
        }

// ----------------------------------------------------------------------------

        /**
         * Private function calculates the mMySubdomainOffset lookup table
         *
         * @return void
         */
        void
        calculate_subdomain_offset()
        {
            Mat< luint > tIJK
                = mBackgroundMesh->get_subdomain_offset_of_proc();

            for( uint l=0; l<gMaxNumberOfLevels; ++l )
            {
                for( uint k=0; k<N; ++k )
                {
                    mMySubdomainOffset[ l ][ k ] = P*tIJK( k, l );
                }
            }
        }

// ----------------------------------------------------------------------------

         /**
          * calculates XZY coordinates for each node
          *
          * @return void
          */
         void
         calculate_node_coordinates()
         {
             // get domain dimensions from settings
             Mat< real > tDomainDimensions
                 = mParameters->get_domain_dimensions();

             // get number of elements on coarsest level from settings
             Mat< luint > tNumberOfElements
                 = mParameters->get_number_of_elements_per_dimension();

             // calculate step width
             real tDeltaX[ gMaxNumberOfLevels ][ N ];

             // calculate width for first level
             for( uint k=0; k<N; ++k )
             {
                 tDeltaX[ 0 ][ k ] = tDomainDimensions( k ) /
                         ( ( real ) ( P * tNumberOfElements( k ) ) );
             }

             // loop over all higher levels
             for( uint l=1; l<gMaxNumberOfLevels; ++l )
             {
                 for( uint k=0; k<N; ++k )
                 {
                     tDeltaX[ l ][ k ] = 0.5*tDeltaX[ l-1 ][ k ];
                 }
             }

             // get domain offset
             Mat< real > tParametersOffset
                 = mParameters->get_domain_offset();

             // domain offset
             real tOffset[ N ];

             // get coords from background mesh
             Mat<real> tOffsetCoords = mBackgroundMesh->get_domain_offset();

             // unflatten coords to a normal array
             for( uint k=0; k<N; ++k )
             {
                 tOffset[ k ] = tOffsetCoords( k );
             }

             // loop over all nodes
             for( auto tNode : mAllBasisOnProc )
             {
                 // get ijk position of node
                 const luint* tIJK = tNode->get_ijk();

                 // get level of node
                 luint tLevel = tNode->get_level();

                 // array containing coordinate
                 real tXYZ[ N ];

                 // loop over all dimensions
                 for( uint k=0; k<N; ++k )
                 {
                     tXYZ[ k ] =
                             ( ( real ) ( tIJK[ k ]
                              + mMySubdomainOffset[ tLevel ][ k ] ) )
                             * tDeltaX[ tLevel ][ k ] + tOffset[ k ];

                 }

                 // write XYZ coordinate into node
                 tNode->set_xyz( tXYZ );
             }
         }

// ----------------------------------------------------------------------------
    };
// ----------------------------------------------------------------------------

        template < uint N, uint P >
        Element *
        Lagrange_Mesh< N, P >::create_element(
                Background_Element_Base* aElement )
        {
            MORIS_ERROR( false, "Don't know how to create Lagrange element.");
            return nullptr;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 2, 1 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 2, 4 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 2, 2 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 2, 9 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 2, 3 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 2, 16 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 2, 4 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 2, 25 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 2, 5 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 2, 36 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 3, 1 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
            = new Lagrange_Element< 3, 8 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 3, 2 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 3, 27 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 3, 3 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 3, 64 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 3, 4 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
                = new Lagrange_Element< 3, 125 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element *
        Lagrange_Mesh< 3, 5 >::create_element(
                Background_Element_Base* aElement )
        {
            Element * aLagrangeElement
            = new Lagrange_Element< 3, 216 >( aElement, mActivePattern );

            return aLagrangeElement;
        }

// ----------------------------------------------------------------------------


    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_ */
