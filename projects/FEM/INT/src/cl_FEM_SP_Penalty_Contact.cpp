/*
 * cl_FEM_SP_Penalty_Contact.cpp
 *
 *  Created on: Feb 13, 2020
 *      Author: ritzert
 */

#include "cl_FEM_SP_Penalty_Contact.hpp"   //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src

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
    SP_Penalty_Contact::SP_Penalty_Contact()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = SP_Property_Type::MATERIAL;
        }


//------------------------------------------------------------------------------
        void SP_Penalty_Contact::reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mMasterVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::MASTER );
            mSlaveVolume  = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::SLAVE );
        }

//------------------------------------------------------------------------------
        void SP_Penalty_Contact::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_Penalty_Contact::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        
//------------------------------------------------------------------------------
        void SP_Penalty_Contact::eval_SP()
        {
            moris::real tEMaster = mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );
            moris::real tESlave = mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );

//std::cout << "tEMaster, tESlave: " << tEMaster << "," << tESlave << std::endl;
//std::cout << "mElementSize" << mElementSize << std::endl;

            mPPVal = mParameters( 0 ) * (mMasterVolume /( mMasterVolume / tEMaster + mSlaveVolume / tESlave)
                    + mSlaveVolume /( mMasterVolume / tEMaster + mSlaveVolume / tESlave)) / mElementSize;
// print(mPPVal,"mPPVal");
//            mPPVal = mParameters( 0 );

        }

////------------------------------------------------------------------------------
//        void SP_Penalty_Contact::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
//        {
//            // get the dof type as a uint
//            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );
//
//            // get the dof type index
//            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );
//
//            // reset the matrix
//            mdPPdMasterDof( tDofIndex ).set_size( 1, mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->get_number_of_space_time_coefficients(), 0.0 );
//
//            // if indirect dependency on the dof type
//            if ( mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->check_dof_dependency( aDofTypes ) )
//            {
//                // compute derivative with indirect dependency through properties
//                mdPPdMasterDof( tDofIndex ).matrix_data()
//                += mParameters( 0 )( 0 )/ (sqrt(mParameters( 1 )( 0 ) * mParameters( 2 )( 0 ) ))
//                    * mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) ) -> dPropdDOF(aDofTypes);
//            }
//        }
//
////------------------------------------------------------------------------------
//        void SP_Penalty_Contact::eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
//        {
//            // get the dof type as a uint
//            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );
//
//            // get the dof type index
//            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );
//
//            // reset the matrix
//            mdPPdMasterDof( tDofIndex ).set_size( 1, mMasterDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
//
//            // if indirect dependency on the dof type
//            if ( mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->check_dof_dependency( aDofTypes ) )
//            {
//                // compute derivative with indirect dependency through properties
//                mdPPdSlaveDof( tDofIndex ).matrix_data()
//                += mParameters( 0 )( 0 )/ (sqrt( mParameters( 1 )( 0 ) * mParameters( 2 )( 0 ) ))
//                    * mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) ) -> dPropdDOF(aDofTypes);
//            }
//        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


