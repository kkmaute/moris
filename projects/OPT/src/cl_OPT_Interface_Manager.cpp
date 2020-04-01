//
// Created by christopherson on 3/20/20.
//

#include "cl_OPT_Interface_Manager.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Interface_Manager::Interface_Manager(ParameterList aParameterList, Cell<std::shared_ptr<Interface>> aInterfaces) : mInterfaces(aInterfaces)
        {
            mNumInterfaces = aInterfaces.size();
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Manager::initialize_advs()
        {
            // Set up ADVs
            Matrix<DDRMat> tGlobalADVs = mInterfaces(0)->initialize_advs();
            uint tCurrentGlobalADVs = tGlobalADVs.length();
            Matrix<DDRMat> tLocalADVs;

            // ADVs per interface
            mNumADVsPerInterface.set_size(mNumInterfaces, 1, 0);
            mNumADVsPerInterface(0) = tCurrentGlobalADVs;

            // ADVs are not shared
            if (!mSharedADVs)
            {
                // Loop through local ADVs
                for (uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local ADVs
                    tLocalADVs = mInterfaces(tInterfaceIndex)->initialize_advs();
                    mNumADVsPerInterface(tInterfaceIndex) = tLocalADVs.length();

                    // Put into the global ADVs
                    tGlobalADVs.resize(tCurrentGlobalADVs + mNumADVsPerInterface(tInterfaceIndex), 1);
                    for (uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface(tInterfaceIndex); tADVIndex++)
                    {
                        tGlobalADVs(tCurrentGlobalADVs + tADVIndex) = tLocalADVs(tADVIndex);
                    }

                    // Update number of global ADVs
                    tCurrentGlobalADVs += mNumADVsPerInterface(tInterfaceIndex);
                }
            }
            return tGlobalADVs;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Manager::get_lower_adv_bounds()
        {
            // Set up ADVs
            Matrix<DDRMat> tGlobalADVs = mInterfaces(0)->get_lower_adv_bounds();
            uint tCurrentGlobalADVs = tGlobalADVs.length();
            Matrix<DDRMat> tLocalADVs;

            // Set upper limit of lower bounds if shared
            if (mSharedADVs)
            {
                for (uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local bounds
                    tLocalADVs = mInterfaces(tInterfaceIndex)->get_lower_adv_bounds();

                    // Compare with current global bounds
                    for (uint tADVIndex = 0; tADVIndex < tCurrentGlobalADVs; tADVIndex++)
                    {
                        tGlobalADVs(tADVIndex) = std::max(tGlobalADVs(tADVIndex), tLocalADVs(tADVIndex));
                    }
                }
            }

            // Loop through local ADVs if not shared
            else
            {
                tGlobalADVs.resize(sum(mNumADVsPerInterface), 1);
                for (uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local ADVs
                    tLocalADVs = mInterfaces(tInterfaceIndex)->get_lower_adv_bounds();

                    // Put into the global ADVs
                    for (uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface(tInterfaceIndex); tADVIndex++)
                    {
                        tGlobalADVs(tCurrentGlobalADVs + tADVIndex) = tLocalADVs(tADVIndex);
                    }

                    // Update number of global ADVs
                    tCurrentGlobalADVs += mNumADVsPerInterface(tInterfaceIndex);
                }
            }
            return tGlobalADVs;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Manager::get_upper_adv_bounds()
        {
            // Set up ADVs
            Matrix<DDRMat> tGlobalADVs = mInterfaces(0)->get_lower_adv_bounds();
            uint tCurrentGlobalADVs = tGlobalADVs.length();
            Matrix<DDRMat> tLocalADVs;

            // Set lower limit of upper bounds if shared
            if (mSharedADVs)
            {
                for (uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local bounds
                    tLocalADVs = mInterfaces(tInterfaceIndex)->get_upper_adv_bounds();

                    // Compare with current global bounds
                    for (uint tADVIndex = 0; tADVIndex < tCurrentGlobalADVs; tADVIndex++)
                    {
                        tGlobalADVs(tADVIndex) = std::min(tGlobalADVs(tADVIndex), tLocalADVs(tADVIndex));
                    }
                }
            }

            // Loop through local ADVs if not shared
            else
            {
                tGlobalADVs.resize(sum(mNumADVsPerInterface), 1);
                for (uint tInterfaceIndex = 1; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local ADVs
                    tLocalADVs = mInterfaces(tInterfaceIndex)->get_upper_adv_bounds();

                    // Put into the global ADVs
                    for (uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface(tInterfaceIndex); tADVIndex++)
                    {
                        tGlobalADVs(tCurrentGlobalADVs + tADVIndex) = tLocalADVs(tADVIndex);
                    }

                    // Update number of global ADVs
                    tCurrentGlobalADVs += mNumADVsPerInterface(tInterfaceIndex);
                }
            }
            return tGlobalADVs;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Manager::get_criteria(Matrix<DDRMat> aNewADVs)
        {
            // Set up global criteria
            Matrix<DDRMat> tGlobalCriteria;
            Matrix<DDRMat> tLocalCriteria;
            uint tCurrentGlobalCriteria = tGlobalCriteria.length();

            // Criteria per interface
            mNumCriteriaPerInterface.set_size(mNumInterfaces, 1);
            mNumCriteriaPerInterface(0) = tCurrentGlobalCriteria;

            // If ADVs are shared, begin new analysis with all ADVs
            if (mSharedADVs)
            {
                for (uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    // Get the local criteria
                    tLocalCriteria = mInterfaces(tInterfaceIndex)->get_criteria(aNewADVs);
                    mNumCriteriaPerInterface(tInterfaceIndex) = tLocalCriteria.length();

                    // Put into the global criteria
                    tGlobalCriteria.resize(tCurrentGlobalCriteria + mNumCriteriaPerInterface(tInterfaceIndex), 1);
                    for (uint tCriteriaIndex = 0; tCriteriaIndex < mNumCriteriaPerInterface(tInterfaceIndex); tCriteriaIndex++)
                    {
                        tGlobalCriteria(tCurrentGlobalCriteria + tCriteriaIndex) = tLocalCriteria(tCriteriaIndex);
                    }

                    // Update number of global criteria
                    tCurrentGlobalCriteria += mNumCriteriaPerInterface(tInterfaceIndex);
                }
            }

            // Otherwise, split the ADVs based on how many are needed per interface
            else
            {
                Matrix<DDRMat> tLocalADVs(0, 0);
                uint tGlobalADVIndex = 0;
                for (uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
                {
                    tLocalADVs.set_size(mNumADVsPerInterface(tInterfaceIndex), 1);
                    for (uint tADVIndex = 0; tADVIndex < mNumADVsPerInterface(tInterfaceIndex); tADVIndex++)
                    {
                        tLocalADVs(tADVIndex) = aNewADVs(tGlobalADVIndex++);
                    }
                    // Get the local criteria
                    tLocalCriteria = mInterfaces(tInterfaceIndex)->get_criteria(tLocalADVs);
                    mNumCriteriaPerInterface(tInterfaceIndex) = tLocalCriteria.length();

                    // Put into the global criteria
                    tGlobalCriteria.resize(tCurrentGlobalCriteria + mNumCriteriaPerInterface(tInterfaceIndex), 1);
                    for (uint tCriteriaIndex = 0; tCriteriaIndex < mNumCriteriaPerInterface(tInterfaceIndex); tCriteriaIndex++)
                    {
                        tGlobalCriteria(tCurrentGlobalCriteria + tCriteriaIndex) = tLocalCriteria(tCriteriaIndex);
                    }

                    // Update number of global criteria
                    tCurrentGlobalCriteria += mNumCriteriaPerInterface(tInterfaceIndex);
                }
            }

            return tGlobalCriteria;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Manager::get_dcriteria_dadv()
        {
            // Set up global criteria gradients
            Matrix<DDRMat> tGlobalCriteriaGradients(sum(mNumCriteriaPerInterface), sum(mNumADVsPerInterface), 0.0);
            Matrix<DDRMat> tLocalCriteriaGradients;
            uint tCurrentGlobalCriteria = 0;
            uint tCurrentGlobalADVs = 0;
            uint tCurrentLocalADVs = 0;

            // Get local criteria gradients from each interface and append accordingly
            for (uint tInterfaceIndex = 0; tInterfaceIndex < mNumInterfaces; tInterfaceIndex++)
            {
                // Get the local criteria
                tLocalCriteriaGradients = mInterfaces(tInterfaceIndex)->get_dcriteria_dadv();
                if (mSharedADVs)
                {
                    tCurrentLocalADVs = tLocalCriteriaGradients.n_cols();
                }
                else
                {
                    tCurrentLocalADVs = mNumADVsPerInterface(tInterfaceIndex);
                }

                // Put into the global criteria
                tGlobalCriteriaGradients({tCurrentGlobalCriteria, tCurrentGlobalCriteria + mNumCriteriaPerInterface(tInterfaceIndex) - 1},
                        {tCurrentGlobalADVs, tCurrentGlobalADVs + tCurrentLocalADVs - 1})
                         = tLocalCriteriaGradients({0, mNumCriteriaPerInterface(tInterfaceIndex) - 1},
                                 {0, tCurrentLocalADVs - 1});

                // Update number of global critera/advs
                tCurrentGlobalCriteria += mNumCriteriaPerInterface(tInterfaceIndex);
                if (!mSharedADVs)
                {
                    tCurrentGlobalADVs += mNumADVsPerInterface(tInterfaceIndex);
                }
            }
            return tGlobalCriteriaGradients;
        }

        // -------------------------------------------------------------------------------------------------------------

    }   // namespace opt
}   // namespace moris