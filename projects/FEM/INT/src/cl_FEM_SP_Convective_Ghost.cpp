
//FEM/INT/src
#include "cl_FEM_SP_Convective_Ghost.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        SP_Convective_Ghost::SP_Convective_Ghost()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ] = Property_Type::DENSITY;
        }

        //------------------------------------------------------------------------------
        void SP_Convective_Ghost::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_Convective_Ghost::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Convective_Ghost::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Convective_Ghost::set_dof_type_list - unknown master slave type." );
                    break;
            }
        }

        //------------------------------------------------------------------------------
        void SP_Convective_Ghost::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "SP_Convective_Ghost::set_property - Unknown aPropertyString : ") +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end() , tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Convective_Ghost::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the density property
            std::shared_ptr< Property > tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get absolute value of u.n
            //            Matrix< DDRMat > tAbs = trans( tVelocityFI->val() ) * mNormal;
            //            real tAbsReal = std::abs( tAbs( 0 ) );
            real tEpsilon = 1e-12;
            Matrix< DDRMat > tNormalDispl = trans( tVelocityFI->val() ) * mNormal;
            real tAbsReal = std::sqrt( tNormalDispl( 0 ) * tNormalDispl( 0 ) + std::pow( tEpsilon, 2.0 ) ) - tEpsilon;

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tDensityProp->val()( 0 ) * tAbsReal * std::pow( mElementSize, 2.0 );
        }

        //------------------------------------------------------------------------------
        void SP_Convective_Ghost::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size(
                    1,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get absolute value of u.n
            //            Matrix< DDRMat > tAbs = trans( tVelocityFI->val() ) * mNormal;
            //            real tAbsReal = std::abs( tAbs( 0 ) );
            real tEpsilon = 1e-12;
            Matrix< DDRMat > tNormalDispl = trans( tVelocityFI->val() ) * mNormal;
            real tAbsReal = std::sqrt( tNormalDispl( 0 ) * tNormalDispl( 0 ) + std::pow( tEpsilon, 2.0 ) ) - tEpsilon;

            // get the density property
            std::shared_ptr< Property > tDensityProp =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                //                // get the sign of u.n
                //                real tSign = 1.0;
                //                if ( tAbs( 0 ) < 0.0 )
                //                {
                //                    tSign = -1.0;
                //                }
                //
                //                // compute contribution from velocity
                //                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                //                        mParameters( 0 ) * std::pow( mElementSize, 2.0 ) * tDensityProp->val()( 0 ) *
                //                        tSign * trans( mNormal ) * tVelocityFI->N();

                // compute contribution from velocity
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        mParameters( 0 ) * std::pow( mElementSize, 2.0 ) * tDensityProp->val()( 0 ) *
                        tNormalDispl( 0 ) * trans( mNormal ) * tVelocityFI->N() /
                        ( tAbsReal + tEpsilon );
            }

            // if density depends on dof
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution from density
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        mParameters( 0 ) * std::pow( mElementSize, 2.0 ) * tAbsReal *
                        tDensityProp->dPropdDOF( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


