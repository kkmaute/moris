/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Set_User_Info.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SET_USER_INFO_HPP_
#define SRC_FEM_CL_FEM_SET_USER_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include "cl_FEM_IWG.hpp"                 //FEM/MSI/src
#include "cl_FEM_IQI.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
        class IWG;
        class IQI;

        //------------------------------------------------------------------------------
        /**
         * Set_User_Info
         */
        class Set_User_Info
        {
            protected :

                // cell of IWG pointers
                moris::Cell< std::shared_ptr< IWG > > mIWGs;

                // cell of IQIs pointers
                moris::Cell< std::shared_ptr< IQI > > mIQIs;

                // set mesh index
                uint mMeshIndex;

                // mesh set name
                std::string mMeshSetName;

                // bool for time sideset
                bool mTimeContinuity = false;

                // bool for time boundary
                bool mTimeBoundary = false;

                // bool for sensitivity analysis computation type
                bool mIsAnalyticalSA = false;

                // enum for FD scheme used for FD SA
                fem::FDScheme_Type mFDSchemeForSA = fem::FDScheme_Type::UNDEFINED;

                // real for finite difference perturbation size
                real mFDPerturbation;

                // bool for forward analysis computation type
                bool mIsAnalyticalFA = true;

                // enum for FD scheme used for FD for forward analysis
                fem::FDScheme_Type mFDSchemeForFA = fem::FDScheme_Type::UNDEFINED;

                // real for finite difference perturbation size for forward analysis
                real mFDPerturbationFA;

               // enum for perturbation strategy used for FD (FA and SA)
                fem::Perturbation_Type mPerturbationStrategy = fem::Perturbation_Type::RELATIVE;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * trivial constructor
                 */
                Set_User_Info(){};

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~Set_User_Info(){};

                //------------------------------------------------------------------------------
                /**
                 * print names
                 */
                void print_names()
                {
                    // print the mesh set name
                    std::cout<<"Mesh set name: "<<mMeshSetName<<std::endl;

                    // print the bool for time sideset
                    std::cout<<"Bool for time sideset: "<<mTimeContinuity<<std::endl;

                    // print IWG names
                    for ( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                    {
                        std::cout<<"IWG name: "<<mIWGs( iIWG )->get_name()<<std::endl;
                    }

                    // print IQI names
                    for ( uint iIQI = 0; iIQI < mIQIs.size(); iIQI++ )
                    {
                        std::cout<<"IQI name: "<<mIQIs( iIQI )->get_name()<<std::endl;
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * set the mesh index
                 * @param[ in ] aMeshIndex mesh index for set
                 */
                void set_mesh_index( uint aMeshIndex )
                {
                    mMeshIndex = aMeshIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the mesh index
                 * @param[ out ] mMeshIndex mesh index for set
                 */
                uint get_mesh_index()
                {
                    return mMeshIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the mesh set name
                 * @param[ in ] aMeshSetName mesh set name
                 */
                void set_mesh_set_name( std::string aMeshSetName )
                {
                    mMeshSetName = aMeshSetName;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the mesh index
                 * @param[ out ] mMeshSetName mesh set name
                 */
                std::string get_mesh_set_name()
                {
                    return mMeshSetName;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the time continuity bool
                 * @param[ in ] aTimeContinuity bool for time sideset
                 */
                void set_time_continuity( bool aTimeContinuity )
                {
                    mTimeContinuity = aTimeContinuity;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the time continuity bool
                 * @param[ out ] mTimeContinuity bool for time sideset
                 */
                const
                bool & get_time_continuity() const
                {
                    return mTimeContinuity;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the time boundary bool
                 * @param[ in ] aTimeBoundary bool for time boundary integral
                 */
                void set_time_boundary( bool aTimeBoundary )
                {
                    mTimeBoundary = aTimeBoundary;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the time boundary bool
                 * @param[ out ] mTimeBoundary bool for time boundary integral
                 */
                const
                bool & get_time_boundary() const
                {
                    return mTimeBoundary;
                }

                //------------------------------------------------------------------------------
                /**
                 * set flag for sensitivity analysis on the set (analytical or finite difference)
                 * @param[ in ] aIsAnalyticalSA bool for sensitivity analysis computation type
                 */
                void set_is_analytical_sensitivity_analysis( bool aIsAnalyticalSA )
                {
                    mIsAnalyticalSA = aIsAnalyticalSA;
                }

                //------------------------------------------------------------------------------
                /**
                 * get flag for sensitivity analysis on the set (analytical or finite difference)
                 * @param[ in ] aIsAnalyticalSA bool for sensitivity analysis computation type
                 */
                const
                bool & get_is_analytical_sensitivity_analysis() const
                {
                    return mIsAnalyticalSA;
                }

                //------------------------------------------------------------------------------
                /**
                 * set FD scheme enum for sensitivity analysis on the set
                 * @param[ in ] aFDSchemeForSA enum for FD scheme used for
                 */
                void set_finite_difference_scheme_for_sensitivity_analysis(
                        enum fem::FDScheme_Type aFDSchemeForSA )
                {
                    mFDSchemeForSA = aFDSchemeForSA;
                }

                //------------------------------------------------------------------------------
                /**
                 * get enum for FD scheme for sensitivity analysis on the set
                 * @param[ out ] mFDSchemeForSA enum for FD scheme used for
                 */
                enum fem::FDScheme_Type get_finite_difference_scheme_for_sensitivity_analysis() const
                {
                    return mFDSchemeForSA;
                }

                //------------------------------------------------------------------------------
                /**
                 * set perturbation size for finite difference
                 * @param[ in ] aFDPerturbation perturbation size
                 */
                void set_finite_difference_perturbation_size( real aFDPerturbation )
                {
                    mFDPerturbation = aFDPerturbation;
                }

                //------------------------------------------------------------------------------
                /**
                 * get perturbation size for finite difference
                 * @param[ out ] mFDPerturbation perturbation size
                 */
                const
                real & get_finite_difference_perturbation_size() const
                {
                    return mFDPerturbation;
                }

                //------------------------------------------------------------------------------
                /**
                 * set flag for forward analysis on the set (analytical or finite difference)
                 * @param[ in ] aIsAnalyticalFA bool for forward analysis computation type
                 */
                void set_is_analytical_forward_analysis( bool aIsAnalyticalFA )
                {
                    mIsAnalyticalFA = aIsAnalyticalFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * get flag for forward analysis on the set (analytical or finite difference)
                 * @param[ in ] aIsAnalyticalFA bool for forward analysis computation type
                 */
                const
                bool & get_is_analytical_forward_analysis() const
                {
                    return mIsAnalyticalFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * set FD scheme enum for forward analysis on the set
                 * @param[ in ] aFDSchemeForFA enum for FD scheme used for
                 */
                void set_finite_difference_scheme_for_forward_analysis(
                        enum fem::FDScheme_Type aFDSchemeForFA )
                {
                    mFDSchemeForFA = aFDSchemeForFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * get enum for FD scheme for forward analysis on the set
                 * @param[ out ] mFDSchemeForFA enum for FD scheme used for
                 */
                enum fem::FDScheme_Type get_finite_difference_scheme_for_forward_analysis() const
                {
                    return mFDSchemeForFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * set perturbation size for finite difference for forward analysis
                 * @param[ in ] aFDPerturbationFA perturbation size
                 */
                void set_finite_difference_perturbation_size_for_forward_analysis( real aFDPerturbationFA )
                {
                    mFDPerturbationFA = aFDPerturbationFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * get perturbation size for finite difference for forward analysis
                 * @param[ out ] mFDPerturbationFA perturbation size
                 */
                const
                real & get_finite_difference_perturbation_size_for_forward_analysis() const
                {
                    return mFDPerturbationFA;
                }

                //------------------------------------------------------------------------------
                /**
                 * set perturbation strategy enum for sensitivity analysis on the set
                 * @param[ in ] aPerturbationStrategy enum for perturbation strategy used for
                 * FA and SA if computed by FD
                 */
                void set_perturbation_strategy(
                        enum fem::Perturbation_Type aPerturbationStrategy )
                {
                    mPerturbationStrategy = aPerturbationStrategy;
                }

                //------------------------------------------------------------------------------
                /**
                 * get perturbation strategy enum for sensitivity analysis on the set
                 * @param[ out ] aPerturbationStrategy enum for perturbation strategy used for
                 * FA and SA if computed by FD
                 */
                enum fem::Perturbation_Type get_perturbation_strategy() const
                {
                    return mPerturbationStrategy;
                }

                //------------------------------------------------------------------------------
                /**
                 * set IWGs
                 * @param[ in ] aIWGs list of IWG pointers
                 */
                void set_IWGs( const moris::Cell< std::shared_ptr< fem::IWG > > & aIWGs )
                {
                    mIWGs = aIWGs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set IWG
                 * @param[ in ] aIWG IWG pointer
                 */
                void set_IWG( std::shared_ptr< fem::IWG > aIWG )
                {
                    mIWGs.push_back( aIWG );
                }

                //------------------------------------------------------------------------------
                /**
                 * get IWGs
                 * @param[ out ] mIWGs list of IWG pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::IWG > > & get_IWGs() const
                {
                    return mIWGs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set IQIs
                 * @param[ in ] aIQIs list of IQI pointers
                 */
                void set_IQIs( const moris::Cell< std::shared_ptr< fem::IQI > > & aIQIs )
                {
                    mIQIs = aIQIs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set IQI
                 * @param[ in ] aIQI IQI pointer
                 */
                void set_IQI( std::shared_ptr< fem::IQI > aIQI )
                {
                    mIQIs.push_back( aIQI );
                }

                //------------------------------------------------------------------------------
                /**
                 * get IQIs
                 * @param[ out ] mIQIs list of IQI pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::IQI > > & get_IQIs() const
                {
                    return mIQIs;
                }

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_SET_USER_INFO_HPP_ */

