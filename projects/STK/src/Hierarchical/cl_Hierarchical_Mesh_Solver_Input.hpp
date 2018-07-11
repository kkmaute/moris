/*
 * cl_Hierachical_Mesh_Solver_Input.hpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_SOLVER_INPUT_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_SOLVER_INPUT_HPP_

#include "linalg.hpp"

#include "cl_Solver_Input.hpp" // DLA/src
#include "cl_Hierarchical_Mesh_Main.hpp"

namespace moris
{
//class Hierarchical_Mesh_Main;

class Hierarchical_Mesh_Solver_Input : public Solver_Input
{
	Hierarchical_Mesh_Main*   mHMR;

    //FIXME get rid of mMat
    moris::Mat< int >         mMat;
    moris::Mat< moris::real > mMatm;

    moris::real mRhsScaling;
    moris::real mRhsOffset;

public :
    // ----------------------------------------------------------------------------------------------
    Hierarchical_Mesh_Solver_Input( moris::real  aRhsScaling,
                                    moris::real  aRhsOffset,
                                    Hierarchical_Mesh_Main*   aHMR );

    // ----------------------------------------------------------------------------------------------
    ~Hierarchical_Mesh_Solver_Input();

    // ----------------------------------------------------------------------------------------------
    // local dimension of the problem
    uint get_num_my_dofs();

    // ----------------------------------------------------------------------------------------------
    // local-to-global map
    Mat <int> get_my_local_global_map();

    // ----------------------------------------------------------------------------------------------
    // element dofs
    uint get_num_element_dof();

    // ----------------------------------------------------------------------------------------------
    // number of elements on proc
    uint get_num_my_elements();

    // ----------------------------------------------------------------------------------------------

//    moris::Mat< moris::real > const & get_element_matrix()
//    {
//        moris::Mat< moris::real > mMatm;
//        return mMatm;
//    };

    // ----------------------------------------------------------------------------------------------
    void get_element_matrix(const uint  & aMyElementInd,
                            Mat< real > & aElementMatrix);

    // ----------------------------------------------------------------------------------------------
//     moris::Mat< int >         const &  get_element_topology()
//   {
//       moris::Mat< int > mMat;
//       return mMat;
//   };

    // ----------------------------------------------------------------------------------------------
    void  get_element_topology(const uint & aMyElementInd,
                               Mat< int > & aElementTopology);

    // ----------------------------------------------------------------------------------------------
    Mat< uint > get_constr_dof();

    // ----------------------------------------------------------------------------------------------

//    moris::Mat< moris::real >          get_element_rhs()
//   {
//       moris::Mat< moris::real > mMatm;
//       return mMatm;
//   };
    // ----------------------------------------------------------------------------------------------
    void get_element_rhs(const uint            & aMyElementInd,
                         Mat< real >           & aElementRHS );
};
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_SOLVER_INPUT_HPP_ */
