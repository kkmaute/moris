#ifndef MORIS_OPTIMIZATION_CL_OPT_ALG_GCMMA_HPP_
#define MORIS_OPTIMIZATION_CL_OPT_ALG_GCMMA_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_OPT_Algorithm.hpp" // Base class // OPT/src

class OptAlgGCMMA : public moris::opt::Algorithm
{

    private:
        moris::uint mResFlag; // Flag from GCMMA describing result of optimization algorithm

        /**
         * @brief External function call for computing objective and constraints, to
         *        interface with gcmma library
         */
        friend void opt_alg_gcmma_func_wrap(
                OptAlgGCMMA* aOptAlgGCMMA,
                int          aIter,
                double*      aAdv,
                double&      aObjval,
                double*      aConval );

        /**
         * @brief External function call for computing sensitivities, to interface
         *        with gcmma library
         */
        friend void opt_alg_gcmma_grad_wrap(
                OptAlgGCMMA* aOptAlgGCMMA,
                double*      aAdv,
                double*      aD_Obj,
                double**     aD_Con,
                int*         aActive );

    public:

        /**
         * Constructor
         */
        OptAlgGCMMA();

        /**
         * Destructor
         */
        ~OptAlgGCMMA();

        /**
         * @brief copy constructor through cloning
         */
        moris::opt::Algorithm*
        clone() const
        {
            return new OptAlgGCMMA(*this );
        }

        /**
         * @brief MORIS interface for solving of optimization problem using
         *        GCMMA
         *
         * @param[in] aOptProb Object of type Problem containing relevant
         *            data regarding ADVs, the objective and constraints
         */
        void solve( moris::opt::Problem* aOptProb );

        /**
         *@brief Prints result of the GCMMA algorithm based on mStopFlag
         */
        void printresult( );
};

#endif /* MORIS_OPTIMIZATION_CL_OPT_ALG_GCMMA_HPP_ */
