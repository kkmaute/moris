/*
 * cl_Solver_Input.hpp
 *
 *  Created on: Apr 6, 2018
 *      Author: schmidt
 */

#ifndef SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_
#define SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
class Dist_Vector;
class Solver_Input
{
public:
    /** Destructor */
    virtual ~Solver_Input(){};

    virtual void set_solution_vector( Dist_Vector * aSolutionVector ){ MORIS_ERROR( false, "Solver_Input::set_solution_vector: not set."); };

    // local dimension of the problem
    virtual moris::uint               get_num_my_dofs()         =0;
    // local-to-global map
    virtual moris::Matrix< DDSMat >          get_my_local_global_map() =0;

    virtual moris::Matrix< DDSMat >          get_my_local_global_overlapping_map( )
    {
        moris::Matrix< DDSMat > aMat;
        //MORIS_ERROR( false, "Solver_Input::get_my_local_global_overlapping_map(): Virtual class not overwritten" );
        return aMat;
    };

    // element dofs
    virtual moris::uint               get_num_element_dof()     =0;
    // number local elements
    virtual moris::uint               get_num_my_elements()     =0;

    virtual moris::Matrix< DDUMat > get_constr_dof()          =0;

    // dense element matrix entries
    //virtual moris::Mat< moris::real >  const & get_element_matrix()   =0;
    //virtual moris::Matrix< DDUMat >          const & get_element_topology() =0;
    //virtual moris::Mat< moris::real >          get_element_rhs()      =0;

    virtual void get_element_matrix(const moris::uint             & aMyElementInd,
                                          moris::Matrix< DDRMat > & aElementMatrix) =0;

    virtual void get_element_topology(const moris::uint             & aMyElementInd,
                                            moris::Matrix< DDSMat > & aElementTopology) =0;

    virtual void get_element_rhs(const moris::uint             & aMyElementInd,
                                       moris::Matrix< DDRMat > & aElementRHS) =0;

    virtual void use_matrix_market_files( ) { MORIS_ERROR(false,"error in use_matrix_market_files"); }

    virtual const char* get_matrix_market_path( )
    {
        //assert(0);
        return NULL;
    }
};
}


#endif /* SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_ */
