///*
// * cl_FEM_IQI_J_Integral.cpp
// *
// *  Created on: Feb 11, 2020
// *      Author: sonne
// */
//#include "cl_FEM_Set.hpp"
//#include "cl_FEM_Field_Interpolator_Manager.hpp"
//#include "cl_FEM_IQI_J_Integral.hpp"
//
//namespace moris
//{
//    namespace fem
//    {
////------------------------------------------------------------------------------
//            IQI_J_Integral::IQI_J_Integral()
//            {
//                // set IQI type
//                mIQIType = vis::Output_Type::J_INTEGRAL;
//
//                // set size for the constitutive model pointer cell
//                mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );
//
//                // populate the constitutive map
//                mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;
//            }
////------------------------------------------------------------------------------
//            void IQI_J_Integral::compute_QI( Matrix< DDRMat > & aQI )
//            {
//                // get indices for properties, CM and SP
//                uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
//
//                // evaluate the QI
//                //aQI =
//
//                         // strain energy-density -> ( trans( mMasterCM( tElastLinIsoIndex )->flux() ) * mMasterCM( tElastLinIsoIndex )->strain() );
//            }
////------------------------------------------------------------------------------
//    }   // end fem namespace
//}       // end moris namespace
//
