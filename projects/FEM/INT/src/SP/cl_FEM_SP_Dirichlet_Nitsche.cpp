/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Dirichlet_Nitsche.cpp
 *
 */

#include "cl_FEM_SP_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Dirichlet_Nitsche::SP_Dirichlet_Nitsche()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        SP_Dirichlet_Nitsche::set_parameters( moris::Cell< Matrix< DDRMat > > aParameters )
        {
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize <= 2,
                    "SP_Dirichlet_Nitsche::set_parameters - either no or 1 constant parameter need to be set." );

            // if geometry measure type is defined
            if ( tParamSize == 2 )
            {
                mGeometryFormulation = mParameters( 1 )( 0 );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Master_Slave > >
        SP_Dirichlet_Nitsche::get_cluster_measure_tuple_list()
        {
            if ( mGeometryFormulation == 1 )
            {
                return { mMasterVolumeTuple, mInterfaceSurfaceTuple };
            }

            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Dirichlet_Nitsche::eval_SP()
        {
            real tElementSize = 0;

            if ( mGeometryFormulation == 1 )
            {
                // get master volume cluster measure value
                const real tMasterVolume =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mMasterVolumeTuple ),
                                        std::get< 1 >( mMasterVolumeTuple ),
                                        std::get< 2 >( mMasterVolumeTuple ) )
                                ->val()( 0 );

                // get interface surface cluster measure value
                const real tInterfaceSurface =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mInterfaceSurfaceTuple ),
                                        std::get< 1 >( mInterfaceSurfaceTuple ),
                                        std::get< 2 >( mInterfaceSurfaceTuple ) )
                                ->val()( 0 );

                tElementSize = tMasterVolume / tInterfaceSurface;
            }
            else
            {
                // get element size cluster measure value
                tElementSize =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mElementSizeTuple ),
                                        std::get< 1 >( mElementSizeTuple ),
                                        std::get< 2 >( mElementSizeTuple ) )
                                ->val()( 0 );
            }

            // get the material property
            const std::shared_ptr< Property >& tPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tPropMaterial->val() / tElementSize;
        }

        //------------------------------------------------------------------------------

        void
        SP_Dirichlet_Nitsche::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            real tElementSize = 0;

            if ( mGeometryFormulation == 1 )
            {
                // get master volume cluster measure value
                const real tMasterVolume =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mMasterVolumeTuple ),
                                        std::get< 1 >( mMasterVolumeTuple ),
                                        std::get< 2 >( mMasterVolumeTuple ) )
                                ->val()( 0 );

                // get interface surface cluster measure value
                const real tInterfaceSurface =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mInterfaceSurfaceTuple ),
                                        std::get< 1 >( mInterfaceSurfaceTuple ),
                                        std::get< 2 >( mInterfaceSurfaceTuple ) )
                                ->val()( 0 );

                tElementSize = tMasterVolume / tInterfaceSurface;
            }
            else
            {
                // get element size cluster measure value
                tElementSize =    //
                        mCluster->get_cluster_measure(
                                        std::get< 0 >( mElementSizeTuple ),
                                        std::get< 1 >( mElementSizeTuple ),
                                        std::get< 2 >( mElementSizeTuple ) )
                                ->val()( 0 );
            }

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get FI for derivative dof type
            Field_Interpolator* tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the material property
            const std::shared_ptr< Property >& tPropMaterial =
                    mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // if material property depends on the dof type
            if ( tPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdPPdMasterDof( tDofIndex ) =
                        mParameters( 0 ) * tPropMaterial->dPropdDOF( aDofTypes ) / tElementSize;
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
