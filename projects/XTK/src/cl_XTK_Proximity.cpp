/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Proximity.cpp
 *
 */

#include "cl_XTK_Proximity.hpp"

using namespace moris;

namespace xtk
{

    // ----------------------------------------------------------------------------------

    void
    Proximity::set_value(
            const moris_index              aGeometryIndex,
            const xtk::Geometric_Proximity aProximity )
    {
        // check that the proximity has not already been set as this would violate the policy of single decisions
        MORIS_ERROR(
                mGeomProximity( aGeometryIndex ) == xtk::Geometric_Proximity::UNDEFINED,
                "xtk::Proximity::set_value() - "
                "Proximity wrt. to geometry #%i is already set and cannot be edited anymore.",
                aGeometryIndex );

        // make sure the geometry index is valid
        MORIS_ASSERT(
                mGeomProximity.size() > (uint) aGeometryIndex,
                "xtk::Proximity::set_value() - "
                "Trying to set proximity wrt. geometry #%i but proximity only defined wrt. %lu geometries.",
                aGeometryIndex,
                mGeomProximity.size() );

        // store the proximity
        mGeomProximity( aGeometryIndex ) = aProximity;
    }

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    Proximity::get_value( const moris_index aGeometryIndex ) const
    {
        // check that the proximity has already been set, otherwise this value should not be exposed to the outside
        MORIS_ERROR(
                mGeomProximity( aGeometryIndex ) == xtk::Geometric_Proximity::UNDEFINED,
                "xtk::Proximity::get_value() - "
                "Proximity wrt. to geometry #%i has not been set.",
                aGeometryIndex );

        // make sure the geometry index is valid
        MORIS_ASSERT(
                mGeomProximity.size() > (uint) aGeometryIndex,
                "xtk::Proximity::get_value() - "
                "Trying to set proximity wrt. geometry #%i but proximity only defined wrt. %lu geometries.",
                aGeometryIndex,
                mGeomProximity.size() );

        // return the proximity value
        return mGeomProximity( aGeometryIndex );
    }

    // ----------------------------------------------------------------------------------

    bool
    Proximity::is_value_set( const moris_index aGeometryIndex ) const
    {
        // make sure the geometry index is valid
        MORIS_ASSERT(
                mGeomProximity.size() > (uint) aGeometryIndex,
                "xtk::Proximity::is_value_set() - "
                "Trying to set proximity wrt. geometry #%i but proximity only defined wrt. %lu geometries.",
                aGeometryIndex,
                mGeomProximity.size() );

        // return whether the proximity is set or not
        return ( mGeomProximity( aGeometryIndex ) != xtk::Geometric_Proximity::UNDEFINED );
    }

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    bool
    do_proximity_values_contradict(
            xtk::Geometric_Proximity aProximityValue1,
            xtk::Geometric_Proximity aProximityValue2 )
    {
        // never process undefined proximity values
        MORIS_ASSERT(
                ( aProximityValue1 != xtk::Geometric_Proximity::UNDEFINED ) && ( aProximityValue2 != xtk::Geometric_Proximity::UNDEFINED ),
                "xtk::do_proximity_values_contradict() - At least one proximity value is 'UNDEFINED'." );

        // values contradict if one indicates that inside, while the other indicates outside
        bool tContradiction1 = ( aProximityValue1 == xtk::Geometric_Proximity::INSIDE ) && ( aProximityValue2 == xtk::Geometric_Proximity::OUTSIDE );
        bool tContradiction2 = ( aProximityValue1 == xtk::Geometric_Proximity::OUTSIDE ) && ( aProximityValue2 == xtk::Geometric_Proximity::INSIDE );

        // return whether a contradiction occurs
        return ( tContradiction1 || tContradiction2 );
    }

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    add_two_proximity_values(
            xtk::Geometric_Proximity aProximityValue1,
            xtk::Geometric_Proximity aProximityValue2 )
    {
        // never process undefined proximity values
        MORIS_ERROR(
                ( aProximityValue1 != xtk::Geometric_Proximity::UNDEFINED ) && ( aProximityValue2 != xtk::Geometric_Proximity::UNDEFINED ),
                "xtk::add_two_proximity_values() - At least one proximity value is 'UNDEFINED'." );

        // matching proximities deliver trivial answer
        if ( aProximityValue1 == aProximityValue2 )
        {
            return aProximityValue1;
        }

        // if one value is interface, the other value automatically decides whether inside, outside, or also interface
        else if ( aProximityValue1 == xtk::Geometric_Proximity::INTERFACE )
        {
            return aProximityValue2;
        }
        else if ( aProximityValue2 == xtk::Geometric_Proximity::INTERFACE )
        {
            return aProximityValue1;
        }

        // contradicting values -> throw error
        else
        {
            MORIS_ERROR( false, "xtk::add_two_proximity_values() - Contradicting proximity values; should never be added." );
            return xtk::Geometric_Proximity::UNDEFINED;
        }

    }    // end function: xtk::add_two_proximity_values()

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    decide_proximity_from_parent_proximities( Vector< xtk::Geometric_Proximity > const & aParentProximities )
    {
        // get the number of proximities to derive new proximity from
        uint tNumProximities = aParentProximities.size();

        // initialize proximity as neutral (interface)
        xtk::Geometric_Proximity tOutputProximityValue = xtk::Geometric_Proximity::INTERFACE;

        // go through all proximities and make decision
        for ( uint iProximity = 0; iProximity < tNumProximities; iProximity++ )
        {
            // get the current proximity
            xtk::Geometric_Proximity tCurrentProximityValue = aParentProximities( iProximity );

            // check that the current proximity is valid
            MORIS_ASSERT( tCurrentProximityValue != xtk::Geometric_Proximity::UNDEFINED,
                    "xtk::decide_proximity_from_parent_proximities() - UNDEFINED proximity value detected." );

            // add the two proximities together
            tOutputProximityValue = add_two_proximity_values( tOutputProximityValue, tCurrentProximityValue );
        }

        // return proximity value
        return tOutputProximityValue;
    }

    // ----------------------------------------------------------------------------------

    // FIXME: this function needs to go
    xtk::Geometric_Proximity
    proximity_vote( Vector< xtk::Geometric_Proximity > const & aParentProximities )
    {
        // get the number of proximities to derive new proximity from
        uint tNumProximities = aParentProximities.size();

        // initialize voting ballot
        Vector< uint > tVotingBallot( 2, 0 );

        // go through all proximities and make decision
        for ( uint iProximity = 0; iProximity < tNumProximities; iProximity++ )
        {
            // get the current proximity
            xtk::Geometric_Proximity tCurrentProximityValue = aParentProximities( iProximity );

            // count votes 
            if( tCurrentProximityValue == xtk::Geometric_Proximity::OUTSIDE )
            {
                tVotingBallot( 0 )++;
            }
            else if( tCurrentProximityValue == xtk::Geometric_Proximity::INSIDE )
            {
                tVotingBallot( 1 )++;
            }
        }

        // return proximity value associated with the winning vote
        if( tVotingBallot( 0 ) > tVotingBallot( 1 ) )
        {
            return xtk::Geometric_Proximity::OUTSIDE;
        }
        else if( tVotingBallot( 0 ) < tVotingBallot( 1 ) )
        {
            return xtk::Geometric_Proximity::INSIDE;
        }
        else
        {
            return xtk::Geometric_Proximity::INTERFACE;
        }

    } // end function: xtk::proximity_vote()

    // ----------------------------------------------------------------------------------

    xtk::Geometric_Proximity
    decide_proximity_from_parent_proximities(
            Vector< const Proximity* > const & aParentProximities,
            const moris_index                aGeometryIndex )
    {
        // get the number of proximities to derive new proximity from
        uint tNumProximities = aParentProximities.size();

        // initialize proximity as neutral (interface)
        xtk::Geometric_Proximity tOutputProximityValue = xtk::Geometric_Proximity::INTERFACE;

        // go through all proximities and make decision
        for ( uint iProximity = 0; iProximity < tNumProximities; iProximity++ )
        {
            // get temporary pointer to the proximity
            const Proximity* tCurrentProximity = aParentProximities( iProximity );

            // check that the current proximity is valid
            MORIS_ASSERT(
                    tCurrentProximity->is_value_set( aGeometryIndex ),
                    "xtk::decide_proximity_from_parent_proximities() - "
                    "%i-th proximity in list is not set wrt. to geometry #%i.",
                    iProximity,
                    aGeometryIndex );

            // get the proximity value
            xtk::Geometric_Proximity tCurrentProximityValue = tCurrentProximity->get_value( aGeometryIndex );

            // add the two proximities together
            tOutputProximityValue = add_two_proximity_values( tOutputProximityValue, tCurrentProximityValue );
        }

        // return proximity value
        return tOutputProximityValue;
    }

    // ----------------------------------------------------------------------------------

}    // namespace xtk