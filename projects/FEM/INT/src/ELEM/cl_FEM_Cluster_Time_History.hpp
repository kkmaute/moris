/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_CL_FEM_CLUSTER_TIME_HISTORY_HPP
#define MORIS_CL_FEM_CLUSTER_TIME_HISTORY_HPP

#include <map>
#include <string>
#include <utility>
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

namespace moris::fem
{
    /** Persistent, named cluster history. Assembly resets must not clear this
     * store. Forward snapshots are immutable during the adjoint sweep.
     * Coefficients are spatial-basis by component, so interpolation remains
     * differentiable when geometry sensitivity perturbs quadrature locations.
     */
    class Cluster_Time_History
    {
      private:
        struct State
        {
            real             start;
            Matrix< DDRMat > values;
        };
        std::map< std::string, std::map< real, State > >  mCommitted;
        std::map< std::string, std::pair< real, State > > mTrial;
        using Point = std::pair< uint, uint >;
        std::map< std::string, std::map< real, std::map< Point, Matrix< DDRMat > > > > mAdjoint;

      public:
        void clear()
        {
            mCommitted.clear();
            mTrial.clear();
            clear_adjoint();
        }

        void clear_adjoint() { mAdjoint.clear(); }

        void stage( const std::string& aName, real aStart, real aEnd, const Matrix< DDRMat >& aValues )
        {
            MORIS_ERROR( aEnd > aStart && aValues.numel() > 0,
                    "Cluster time history requires a positive slab duration and nonempty values." );
            mTrial[ aName ] = { aEnd, { aStart, aValues } };
        }

        void rollback( const std::string& aName ) { mTrial.erase( aName ); }

        void commit( const std::string& aName )
        {
            auto tTrial = mTrial.find( aName );
            MORIS_ERROR( tTrial != mTrial.end(), "Cluster time history: no staged state." );
            auto& tHistory = mCommitted[ aName ];
            MORIS_ERROR( tHistory.empty() || tHistory.rbegin()->first == tTrial->second.second.start,
                    "Cluster time history: slabs must be committed once in chronological order." );
            tHistory.emplace( tTrial->second.first, std::move( tTrial->second.second ) );
            mTrial.erase( tTrial );
        }

        const Matrix< DDRMat >& values( const std::string& aName, real aEnd ) const
        {
            auto tVariable = mCommitted.find( aName );
            MORIS_ERROR( tVariable != mCommitted.end(), "Cluster time history: unknown variable." );
            auto tState = tVariable->second.find( aEnd );
            MORIS_ERROR( tState != tVariable->second.end(), "Cluster time history: missing slab." );
            return tState->second.values;
        }

        real start( const std::string& aName, real aEnd ) const
        {
            (void)values( aName, aEnd );
            return mCommitted.at( aName ).at( aEnd ).start;
        }

        void set_adjoint( const std::string& aName, real aTime, Point aPoint, const Matrix< DDRMat >& aValue )
        {
            mAdjoint[ aName ][ aTime ][ aPoint ] = aValue;
        }

        Matrix< DDRMat > take_adjoint( const std::string& aName, real aTime, Point aPoint )
        {
            // A missing contribution is a read-only lookup: reverse assembly
            // may request points/times for which no carry was ever staged.
            auto tVariable = mAdjoint.find( aName );
            if ( tVariable == mAdjoint.end() )
            {
                return {};
            }

            auto tTime = tVariable->second.find( aTime );
            if ( tTime == tVariable->second.end() )
            {
                return {};
            }

            auto tEntry = tTime->second.find( aPoint );
            if ( tEntry == tTime->second.end() )
            {
                return {};
            }

            Matrix< DDRMat > tValue = std::move( tEntry->second );
            tTime->second.erase( tEntry );

            // Prune only consumed adjoint storage. Other points/times and the
            // committed forward trajectory remain available to the reverse sweep.
            if ( tTime->second.empty() )
            {
                tVariable->second.erase( tTime );

                if ( tVariable->second.empty() )
                {
                    mAdjoint.erase( tVariable );
                }
            }

            return tValue;
        }
    };
}    // namespace moris::fem
#endif
