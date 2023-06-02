/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Facet.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_BACKGROUND_FACET_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_BACKGROUND_FACET_HPP_

#include "HMR_Globals.hpp" //HMR/src
#include "typedefs.hpp"
#include "cl_Bitset.hpp" //CNT/src

namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------

        class Background_Element_Base;

//-----------------------------------------------------------------------------
        class Background_Facet
        {
            //! Pointer to parent. Points to null for faces
            //! on coarsest level.
            // Background_Facet *      mParent;

            //! Contains the ID of the proc that owns the element.
            // moris_id                    mOwner;

            //! Tells if a face is active
            // Bitset< gNumberOfPatterns > mActiveFlags;

            //! Tells if a face is refined
            // Bitset< gNumberOfPatterns > mRefinedFlags;

            //! Tells if the face has children.
            // bool                        mChildrenFlag = false;

            //! index in memory
            // uint                        mMemoryIndex;

            //! reference element ( the element with the lower id )
            Background_Element_Base *   mLeaderElement;
            Background_Element_Base *   mFollowerElement;

            //! face index for leader element
            uint                        mIndexOnLeader;

            //! multi purpose flag
            bool                        mFlag = false;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * Default constructor
             */
            Background_Facet(       Background_Element_Base * aElementA,
                                    Background_Element_Base * aElementB,
                              const uint                    & aIndexOnElementA  );

//------------------------------------------------------------------------------

            Background_Facet(       Background_Element_Base * aElement,
                              const uint                    & aIndexOnElement );

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Background_Facet(){};

//------------------------------------------------------------------------------

            /**
             * sets the flag
             */
            void flag();

//------------------------------------------------------------------------------

            /**
             * resets the flag
             */
            void unflag();

//--------------------------------------------------------------------------------

            /**
             * test if flag is set
             */
            bool is_flagged() const;

//--------------------------------------------------------------------------------

            /**
             * get pointer to leader element
             */
            Background_Element_Base * get_leader();
//--------------------------------------------------------------------------------

            /**
             * get pointer to follower element
             */
            Background_Element_Base * get_follower();

//--------------------------------------------------------------------------------

            /**
             * sets the follower element
             */
            void set_follower( Background_Element_Base * aElement );

//--------------------------------------------------------------------------------

            /**
             * returns the index in relation to the leader element
             */
            uint get_index_on_leader() const;

//--------------------------------------------------------------------------------

            /**
             * returns the index in relation to the leader element
             */
            uint get_index_on_follower() const;

//--------------------------------------------------------------------------------

            /**
             * returns the index in relation to the leader element
             */
            uint get_index_on_other( const uint & aIndex ) const;

//--------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_BACKGROUND_FACET_HPP_ */

