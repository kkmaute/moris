//
// Created by christopherson on 2/7/20.
//

#include "cl_OPT_Problem_Rosenbrock.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Problem_Rosenbrock::get_constraint_types()
        {
            Matrix<DDSMat> tConstraintTypes(2, 1, 1);

            return tConstraintTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_objectives()
        {
            Matrix<DDRMat> tObjectives(1, 1);
            tObjectives(0) = (1 - mADVs(0)) * mCriteria(0, 0) + 100 * (mADVs(1) - pow(mADVs(0), 2)) * mCriteria(1, 0);

            return tObjectives;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_constraints()
        {
            Matrix<DDRMat> tConstraints(2, 1, 0.0);
            if (mConstrained)
            {
                tConstraints(0) = pow(mCriteria(0), 2) * (mADVs(0) - 1) - mADVs(1) + 1;
                tConstraints(1) = mADVs(0) + mADVs(1) - 2;
            }

            return tConstraints;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_dobjective_dadv()
        {
            Matrix<DDRMat> tDObjectiveDADV(1, 2);
            tDObjectiveDADV(0) = -mCriteria(0) - 200 * mADVs(0) * mCriteria(1);
            tDObjectiveDADV(1) = 100 * mCriteria(1);

            return tDObjectiveDADV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_dconstraint_dadv()
        {
            Matrix<DDRMat> tDConstraintDADV(2, 2, 0.0);
            if (mConstrained)
            {
                tDConstraintDADV(0, 0) = pow(mCriteria(0), 2);
                tDConstraintDADV(0, 1) = -1;
                tDConstraintDADV(1, 0) = 1;
                tDConstraintDADV(1, 1) = 1;
            }

            return tDConstraintDADV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_dobjective_dcriteria()
        {
            Matrix<DDRMat> tDObjectiveDCriteria(1, 2);
            tDObjectiveDCriteria(0) = 1 - mADVs(0);
            tDObjectiveDCriteria(1) = 100 * (mADVs(1) - pow(mADVs(0), 2));

            return tDObjectiveDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem_Rosenbrock::compute_dconstraint_dcriteria()
        {
            Matrix<DDRMat> tDConstraintDCriteria(2, 2, 0.0);
            if (mConstrained)
            {
                tDConstraintDCriteria(0, 0) = -2 * mCriteria(0) * (mADVs(0) - 1);
                tDConstraintDCriteria(0, 1) = 0;
                tDConstraintDCriteria(1, 0) = 0;
                tDConstraintDCriteria(1, 1) = 0;
            }

            return tDConstraintDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
