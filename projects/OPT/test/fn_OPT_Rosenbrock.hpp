#ifndef MORIS_FN_OPT_ROSENBROCK_HPP
#define MORIS_FN_OPT_ROSENBROCK_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        void initialize_rosenbrock(
                Matrix<DDRMat>& aADVs,
                Matrix<DDRMat>& aLowerBounds,
                Matrix<DDRMat>& aUpperBounds)
        {
            // Initial Guess
            aADVs.set_size(2, 1);
            aADVs(0) = 0.8;
            aADVs(1) = 1.2;

            // Lower Bounds
            aLowerBounds.set_size(2, 1);
            aLowerBounds(0) = -2.0;
            aLowerBounds(1) = -2.0;

            // Upper Bounds
            aUpperBounds.set_size(2, 1);
            aUpperBounds(0) = 2.0;
            aUpperBounds(1) = 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> get_criteria_rosenbrock(const Matrix<DDRMat>& aADVs)
        {
            Matrix<DDRMat> tCriteria(2, 1);
            tCriteria(0) = 1 - aADVs(0);
            tCriteria(1) = aADVs(1) - pow(aADVs(0), 2);

            return tCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> get_dcriteria_dadv_rosenbrock(const Matrix<DDRMat>& aADVs)
        {
            Matrix<DDRMat> tDCriteria(2, 2, 0.0);
            tDCriteria(0, 0) = -1;
            tDCriteria(0, 1) = 0;
            tDCriteria(1, 0) = -2 * aADVs(0);
            tDCriteria(1, 1) = 1;

            return tDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> get_constraint_types_rosenbrock()
        {
            Matrix<DDSMat> tConstraintTypes(2, 1, 1);

            return tConstraintTypes;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_objectives_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tObjectives(1, 1);
            tObjectives(0) = (1 - aADVs(0)) * aCriteria(0) + 100 * (aADVs(1) - pow(aADVs(0), 2)) * aCriteria(1);

            return tObjectives;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_constraints_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tConstraints(2, 1);
            tConstraints(0) = pow(aCriteria(0), 2) * (aADVs(0) - 1) - aADVs(1) + 1;
            tConstraints(1) = aADVs(0) + aADVs(1) - 2;

            return tConstraints;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_dobjective_dadv_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tDObjectiveDADV(1, 2);
            tDObjectiveDADV(0) = -aCriteria(0) - 200 * aADVs(0) * aCriteria(1);
            tDObjectiveDADV(1) = 100 * aCriteria(1);

            return tDObjectiveDADV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_dobjective_dcriteria_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tDObjectiveDCriteria(1, 2);
            tDObjectiveDCriteria(0) = 1 - aADVs(0);
            tDObjectiveDCriteria(1) = 100 * (aADVs(1) - pow(aADVs(0), 2));

            return tDObjectiveDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_dconstraint_dadv_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tDConstraintDADV(2, 2);
            tDConstraintDADV(0, 0) = pow(aCriteria(0), 2);
            tDConstraintDADV(0, 1) = -1;
            tDConstraintDADV(1, 0) = 1;
            tDConstraintDADV(1, 1) = 1;

            return tDConstraintDADV;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> compute_dconstraint_dcriteria_rosenbrock(
                const Matrix<DDRMat>& aADVs,
                const Matrix<DDRMat>& aCriteria)
        {
            Matrix<DDRMat> tDConstraintDCriteria(2, 2);
            tDConstraintDCriteria(0, 0) = -2 * aCriteria(0) * (aADVs(0) - 1);
            tDConstraintDCriteria(0, 1) = 0;
            tDConstraintDCriteria(1, 0) = 0;
            tDConstraintDCriteria(1, 1) = 0;

            return tDConstraintDCriteria;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

#endif //MORIS_FN_OPT_ROSENBROCK_HPP
