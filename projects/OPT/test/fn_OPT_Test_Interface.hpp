/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_Test_Interface.hpp
 *
 */

#ifndef MORIS_FN_OPT_TEST_INTERFACE_HPP
#define MORIS_FN_OPT_TEST_INTERFACE_HPP

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        void initialize_test_1(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            // Initial Guess
            aADVs.set_size(2, 1);
            aADVs(0) = 1.0;
            aADVs(1) = 2.0;

            // Lower Bounds
            aLowerBounds.set_size(2, 1);
            aLowerBounds(0) = 1.0;
            aLowerBounds(1) = 2.0;

            // Upper Bounds
            aUpperBounds.set_size(2, 1);
            aUpperBounds(0) = 1.0;
            aUpperBounds(1) = 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void initialize_test_2(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            // Initial Guess
            aADVs.set_size(2, 1);
            aADVs(0) = 3.0;
            aADVs(1) = 4.0;

            // Lower Bounds
            aLowerBounds.set_size(2, 1);
            aLowerBounds(0) = 3.0;
            aLowerBounds(1) = 4.0;

            // Upper Bounds
            aUpperBounds.set_size(2, 1);
            aUpperBounds(0) = 3.0;
            aUpperBounds(1) = 4.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void initialize_test_3(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            // Initial Guess
            aADVs.set_size(2, 1);
            aADVs(0) = 5.0;
            aADVs(1) = 6.0;

            // Lower Bounds
            aLowerBounds.set_size(2, 1);
            aLowerBounds(0) = 5.0;
            aLowerBounds(1) = 6.0;

            // Upper Bounds
            aUpperBounds.set_size(2, 1);
            aUpperBounds(0) = 5.0;
            aUpperBounds(1) = 6.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void initialize_test_4(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            // Initial Guess
            aADVs.set_size(2, 1);
            aADVs(0) = 7.0;
            aADVs(1) = 8.0;

            // Lower Bounds
            aLowerBounds.set_size(2, 1);
            aLowerBounds(0) = 7.0;
            aLowerBounds(1) = 8.0;

            // Upper Bounds
            aUpperBounds.set_size(2, 1);
            aUpperBounds(0) = 7.0;
            aUpperBounds(1) = 8.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> get_criteria_test(const Matrix<DDRMat>& aADVs)
        {
            Matrix<DDRMat> tCriteria(2, 1);
            tCriteria(0) = aADVs(0);
            tCriteria(1) = aADVs(1);

            return tCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> get_dcriteria_dadv_test(const Matrix<DDRMat>& aADVs)
        {
            Matrix<DDRMat> tDCriteria(2, 2, 0.0);
            tDCriteria(0, 0) = aADVs(0);
            tDCriteria(0, 1) = 0;
            tDCriteria(1, 0) = 0;
            tDCriteria(1, 1) = aADVs(1);

            return tDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

#endif //MORIS_FN_OPT_TEST_INTERFACE_HPP

