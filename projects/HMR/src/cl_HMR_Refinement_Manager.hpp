/*
 * cl_HMR_Refinement_Manager.hpp
 *
 *  Created on: May 30, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_
#define SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_

#include "typedefs.hpp" //COR/src
//#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_Mat.hpp" //LNA/src
//#include "cl_HMR_Parameters.hpp" //HMR/src
//#include "cl_HMR_Parameters.hpp" //HMR/src

//#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
//#include "cl_HMR_Knot.hpp" //HMR/src
//#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src


namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------

        // forward declaration of HMR object
        class Lagrange_Mesh_Base;

//-------------------------------------------------------------------------------
        class Refinement_Manager
        {

            // pointer to mesh object
            Lagrange_Mesh_Base * mMesh;

            //! max surface level
            const uint mMaxSurfaceLevel;

            //! max volume level
            const uint mMaxVolumeLevel;

//-------------------------------------------------------------------------------
    public:
//-------------------------------------------------------------------------------

        /**
         * \brief object that manages refinement
         *
         *  @param[inout]  aMesh     mesh the operation is performed on
         */
        Refinement_Manager ( Lagrange_Mesh_Base * aMesh );

//-------------------------------------------------------------------------------

        /* destructor for refinement manager */
        ~Refinement_Manager(){}

//-------------------------------------------------------------------------------

        void
        flag_against_nodal_field(
                const Mat< real > & aField,
                const real          aLowerBound=0.0,
                const real          aUpperBound=0.0 );

//-------------------------------------------------------------------------------

        void
        flag_against_elemental_field(
                const Mat< real > & aField,
                const real          aLowerBound=0.0);

//-------------------------------------------------------------------------------
        };
//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_ */
