/*
 * cl_FEM_Phase_User_Info.hpp
 *
 *  Created on: Oct 01, 2020
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_
#define SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_

//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
namespace moris
{
    namespace fem
    {

        class Constitutive_Model;
        class Stabilization_Parameter;

        //------------------------------------------------------------------------------
        /**
         * Phase_User_Info
         */
        class Phase_User_Info
        {
            protected :

                // phase name
                std::string mPhaseName;

                // phase index
                sint mPhaseIndex;

                // phase dof type list
                moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

                // phase dv type list
                moris::Cell< moris::Cell< PDV_Type > > mDvTypes;

                // master and slave constitutive models
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mCMs;

                // stabilization parameters
                moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mSPs;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * trivial constructor
                 */
                Phase_User_Info(){};

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~Phase_User_Info(){};

                //------------------------------------------------------------------------------
                /**
                 * print names
                 */
                void print_names()
                {
                    // print the phase name
                    std::cout<<"Phase name: "<<mPhaseName<<std::endl;

                    // print the phase index
                    std::cout<<"Phase index: "<<mPhaseIndex<<std::endl;

                    // print CM names
                    for ( uint iCM = 0; iCM < mCMs.size(); iCM++ )
                    {
                        std::cout<<"CM name: "<<mCMs( iCM )->get_name()<<std::endl;
                    }

                    // print SP names
                    for ( uint iSP = 0; iSP < mSPs.size(); iSP++ )
                    {
                        std::cout<<"SP name: "<<mSPs( iSP )->get_name()<<std::endl;
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * set the phase name
                 * @param[ in ] aPhaseName mesh set name
                 */
                void set_phase_name( std::string aPhaseName )
                {
                    mPhaseName = aPhaseName;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the phase index
                 * @param[ out ] mPhaseName mesh set name
                 */
                std::string get_phase_name()
                {
                    return mPhaseName;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the phase index
                 * @param[ in ] aPhaseIndex phase index
                 */
                void set_phase_index( sint aPhaseIndex )
                {
                    mPhaseIndex = aPhaseIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the phase index
                 * @param[ out ] mPhaseIndex mesh index for set
                 */
                sint get_phase_index()
                {
                    return mPhaseIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * set dof type list
                 * @param[ out ] mDofTypes list of group of dof types
                 */
                void set_dof_type_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
                {
                    mDofTypes = aDofTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * get dof type list
                 * @param[ out ] mDofTypes list of group of dof types
                 */
                const
                moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list()
                {
                    return mDofTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * set dv type list
                 * @param[ out ] mDvTypes list of group of dv types
                 */
                void set_dv_type_list( const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes )
                {
                    mDvTypes = aDvTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * get dv type list
                 * @param[ out ] mDvTypes list of group of dv types
                 */
                const
                moris::Cell< moris::Cell< PDV_Type > > & get_dv_type_list()
                {
                    return mDvTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * set CMs
                 * @param[ in ] aCMs list of CM pointers
                 */
                void set_CMs( const moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & aCMs )
                {
                    mCMs = aCMs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set CM
                 * @param[ in ] aCM constitutive model pointer
                 */
                void set_CM( std::shared_ptr< fem::Constitutive_Model > aCM )
                {
                    mCMs.push_back( aCM );
                }

                //------------------------------------------------------------------------------
                /**
                 * get CMs
                 * @param[ out ] mCMs list of constitutive model pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & get_CMs() const
                {
                    return mCMs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set SPs
                 * @param[ in ] aSPs list of stabilization parameter pointers
                 */
                void set_SPs( std::shared_ptr< fem::Stabilization_Parameter > aSPs )
                {
                    mSPs.push_back( aSPs );
                }

                //------------------------------------------------------------------------------
                /**
                 * get SPs
                 * @param[ out ] mSPs list of stabilization parameter pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & get_SPs() const
                {
                    return mSPs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set SP
                 * @param[ in ] aSP stabilization parameter pointer
                 */
                void set_SP( std::shared_ptr< fem::Stabilization_Parameter > aSP )
                {
                    mSPs.push_back( aSP );
                }

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_ */
