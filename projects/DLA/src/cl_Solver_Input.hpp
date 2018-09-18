/*
 * cl_Solver_Input.hpp
 *
 *  Created on: Apr 6, 2018
 *      Author: schmidt
 */

#ifndef SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_
#define SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_

#include "linalg.hpp"

namespace moris
{
class Solver_Input
{

public:
    /** Destructor */
    virtual ~Solver_Input(){};

    // local dimension of the problem
    virtual moris::uint               get_num_my_dofs()         =0;
    // local-to-global map
    virtual moris::Mat <int>          get_my_local_global_map() =0;

    virtual moris::Mat <int>          get_my_local_global_overlapping_map( )
    {

        moris::Mat< int > aMat;
        //MORIS_ERROR( false, "Solver_Input::get_my_local_global_overlapping_map(): Virtual class not overwritten" );
        return aMat;
    };

    // element dofs
    virtual moris::uint               get_num_element_dof()     =0;
    // number local elements
    virtual moris::uint               get_num_my_elements()     =0;

    virtual moris::Mat< moris::uint > get_constr_dof()          =0;

    // dense element matrix entries
    //virtual moris::Mat< moris::real >  const & get_element_matrix()   =0;
    //virtual moris::Mat< int >          const & get_element_topology() =0;
    //virtual moris::Mat< moris::real >          get_element_rhs()      =0;

    virtual void get_element_matrix(const moris::uint               & aMyElementInd,
                                          moris::Mat< moris::real > & aElementMatrix) =0;

    virtual void get_element_topology(const moris::uint       & aMyElementInd,
                                            moris::Mat< int > & aElementTopology)     =0;

    virtual void get_element_rhs(const moris::uint               & aMyElementInd,
                                       moris::Mat< moris::real > & aElementRHS)       =0;

    virtual void use_matrix_market_files( )  { MORIS_ERROR(false,"error in use_matrix_market_files"); }

    virtual const char* get_matrix_market_path( )
    {
        //assert(0);
        return NULL;
    }
};
}


#endif /* SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_ */
