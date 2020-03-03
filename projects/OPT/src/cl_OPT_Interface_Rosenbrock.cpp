//
// Created by christopherson on 2/7/20.
//

#include "cl_OPT_Interface_Rosenbrock.hpp"

namespace moris
{
    namespace opt
    {
        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Rosenbrock::initialize_advs()
        {
            Matrix<DDRMat> tADVs(2, 1);
            tADVs(0) = 0.8;
            tADVs(1) = 1.2;

            return tADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Rosenbrock::get_lower_adv_bounds()
        {
            Matrix<DDRMat> tLowerBounds(2, 1);
            tLowerBounds(0) = -2.0;
            tLowerBounds(1) = -2.0;

            return tLowerBounds;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Rosenbrock::get_upper_adv_bounds()
        {
            Matrix<DDRMat> tUpperBounds(2, 1);
            tUpperBounds(0) = 2.0;
            tUpperBounds(1) = 2.0;

            return tUpperBounds;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interface_Rosenbrock::begin_new_analysis(Matrix<DDRMat> aNewADVs)
        {
            mADVs = aNewADVs;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Rosenbrock::get_criteria()
        {
            Matrix<DDRMat> tCriteria(2, 1);
            tCriteria(0) = 1 - mADVs(0);
            tCriteria(1) = mADVs(1) - pow(mADVs(0), 2);

            return tCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interface_Rosenbrock::get_dcriteria_dadv()
        {
            Matrix<DDRMat> tDCriteria(2, 2, 0.0);
            tDCriteria(0, 0) = -1;
            tDCriteria(0, 1) = 0;
            tDCriteria(1, 0) = -2 * mADVs(0);
            tDCriteria(1, 1) = 1;

            return tDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

    }   // namespace opt
}   // namespace moris
