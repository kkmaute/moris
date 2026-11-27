/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_CL_FEM_IWG_TIME_DG_HPP
#define MORIS_CL_FEM_IWG_TIME_DG_HPP

#include "cl_FEM_IWG.hpp"

namespace moris::fem
{
    /** Spatial-quadrature inertia/history block of the report's Eq. (54).
     * Register on a time-continuity set. Velocity is stored as spatial
     * coefficients in cluster ISVs and interpolated at each quadrature point.
     */
    class IWG_Time_DG_Mass : public IWG
    {
        enum Property_Index
        {
            DENSITY,
            THICKNESS,
            INITIAL_DISPLACEMENT,
            INITIAL_VELOCITY,
            INITIAL_TIME,
            COUNT
        };
        std::string      history_name() const;
        real             start_time() const;
        real             duration() const;
        bool             initial_slab() const;
        Matrix< DDRMat > initial_coefficients( uint aProperty );
        Matrix< DDRMat > incoming_displacement();
        Matrix< DDRMat > incoming_velocity();
        Matrix< DDRMat > spatial_shape();
        void             assemble( real aWeight, bool aResidual, bool aJacobian, bool aPrevious );

      public:
        IWG_Time_DG_Mass();
        bool has_time_history() const override { return true; }
        void commit_time_history() override;
        void compute_residual( real aWeight ) override;
        void compute_jacobian( real aWeight ) override;
        void compute_jacobian_and_residual( real aWeight ) override;
        void compute_jacobian_previous( real aWeight ) override;
        void compute_dRdp( real aWeight ) override;
        void add_adjoint_history( real aWeight, uint aElement, uint aPoint, const Vector< Matrix< DDRMat > >& aAdjoint, Vector< Matrix< DDRMat > >& aRHS ) override;
    };

    /** Decorates an ordinary conservative spatial IWG with the report's
     * reduced time tests and incoming collocation. Trial interpolation and
     * constitutive tangents are unchanged. Spatial surface and ghost terms
     * receive exactly the same temporal treatment as bulk forces.
     */
    class IWG_Time_DG_Spatial : public IWG
    {
        std::shared_ptr< IWG > mSpatial;
        void                   configure_spatial();
        void                   add_lifted_rows( const Matrix< DDRMat >& aMatrix, uint aTemporalRow, Matrix< DDRMat >& aResult );
        void                   assemble( real aWeight, bool aResidual, bool aJacobian );

      public:
        explicit IWG_Time_DG_Spatial( std::shared_ptr< IWG > aSpatial );
        bool uses_dg_time_tests() const override { return true; }
        void compute_residual( real aWeight ) override;
        void compute_jacobian( real aWeight ) override;
        void compute_jacobian_and_residual( real aWeight ) override;
        void compute_jacobian_previous( real aWeight ) override;
        void compute_dRdp( real aWeight ) override;
    };
}    // namespace moris::fem
#endif
