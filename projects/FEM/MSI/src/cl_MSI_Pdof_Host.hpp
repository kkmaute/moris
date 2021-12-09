/*
 * cl_Pdof_Host.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_PDOF_HOST_HPP_
#define SRC_FEM_CL_PDOF_HOST_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "fn_unique.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MSI_Adof.hpp"

namespace moris
{
    namespace fem
    {
        class Node_Base;
    }
}

namespace moris
{
    namespace MSI
    {
        struct Pdof
        {
                moris::uint        mDofTypeIndex;
                moris::uint        mTimeStepIndex;
                Matrix< DDSMat >   mAdofIds;
                Matrix< DDRMat >   mTmatrix;

                moris::uint        mElementalSolVecEntry;

                moris::Cell < Adof* >  mAdofPtrList;   //FIXME delete this list after call to get adof ids or replace it
        };

        //-------------------------------------------------------------------------------------------------
        class Model_Solver_Interface;
        class Pdof_Host
        {
            private:

                Matrix< DDUMat >                        mPdofTypeExist;         // Vector indicates if dof type exists. FIXME make this a global matrix in dof manager and delete after costruction
                moris::Cell< moris::Cell< Pdof* > >     mListOfPdofTimePerType; // List of all pdofs per time per dof type. outer cell is type, inner cell is time

                void create_adofs_based_on_Tmatrix(
                        const Matrix< DDUMat >               & aTimeLevelOffsets,
                        moris::Cell< moris::Cell< Adof * > > & aAdofListz,
                        Model_Solver_Interface               * aModelSolverInterface);

                void create_adofs_based_on_pdofs(
                        const Matrix< DDUMat >               & aTimeLevelOffsets,
                        moris::Cell< moris::Cell< Adof * > > & aAdofList );

            protected:

                fem::Node_Base * mNodeObj;
                moris_id  mNodeID;
                //FIXME Add interpolation order

            public:

                Pdof_Host() {};

                Pdof_Host(
                        const moris::uint   aNumUsedDofTypes,
                        fem::Node_Base    * aNodeObj );

                ~Pdof_Host();

                /**
                 * @brief Sets a dof type and time to this pdof host. This function is tested by the test [Pdof_host_set_dof_type]
                 *
                 * @param[in] aDof_Type          Dof type
                 * @param[in] aTimeSteps         Time steps.
                 * @param[in] aNumUsedDofTypes   Globally used dof types.
                 * @param[in] aPdofTypeMap       Paw which related a dof type enum to the index.
                 *
                 */
                void set_pdof_type(
                        const enum Dof_Type         aDof_Type,
                        const Matrix< DDUMat >    & aTimePerDofType,
                        const moris::uint           aNumUsedDofTypes,
                        const Matrix< DDSMat >    & aPdofTypeMap);

                /**
                 * @brief Gets the adofs for all the pdofs in this pdof host. This function is tested by the test [Pdof_Host_Get_Adofs]
                 *
                 * @param[in] aTimeLevelOffsets  Offsets for this doftype and time      FIXME
                 * @param[in] aAdofList          List containing all the adofs.
                 *
                 */
                void get_adofs(
                        const Matrix< DDUMat >               & aTimeLevelOffsets,
                        moris::Cell< moris::Cell< Adof * > > & aAdofList,
                        Model_Solver_Interface               * aModelSolverInterface,
                        const bool                           & aUseHMR );

                /**
                 * @brief Gets the adofs Ids for all the pdofs in this pdof host. This function is tested by the test [Pdof_Host_Get_Adofs]
                 *
                 * @param[in] aTimeLevelOffsets       Offsets for this doftype and time
                 * @param[in] aAdofList               List containing all the adofs.
                 * @param[in] aModelSolverInterface   Pointer to the model solver interface
                 * @param[in] aUseHMR                 Bolean which indicates if HMR is used and thus, multiple adofs per pdof.
                 *
                 */
                void get_adofs_ids();

                /**
                 * @brief Set the t-matrix values for all the pdofs. This function is tested by the test [Pdof_Host_Get_Adofs]
                 *
                 */
                void set_t_matrix( const bool & aUseHMR,
                        Model_Solver_Interface * aModelSolverInterface);

                /**
                 * @brief Print t-matrix values for all the pdofs.
                 *
                 * @param[in] aModelSolverInterface   Pointer to the model solver interface
                 * @param[in] aUseHMR                 Bolean which indicates if HMR is used and thus, multiple adofs per pdof.
                 *
                 */
                void print_t_matrix(
                        const bool              & aUseHMR,
                        Model_Solver_Interface * aModelSolverInterface );

                //-------------------------------------------------------------------------------------------------
                // FIXME member function need comments

                fem::Node_Base * get_node_obj_ptr()
                {
                    return mNodeObj;
                }

                void set_pointer_to_Tmatrix()
                {}

                void set_adof_IDs()
                {}

                moris::uint get_num_pdofs();

                moris::Cell< Pdof* > & get_pdof_time_list( const moris::sint & aDofTypeIndex)
                {
                    return mListOfPdofTimePerType( aDofTypeIndex );
                }

                moris::uint get_num_time_levels_of_type( const moris::uint & aDofTypeInd )
                {
                    return mListOfPdofTimePerType( aDofTypeInd ).size();
                }

                moris::Cell< moris::Cell< Pdof* > > & get_pdof_hosts_pdof_list()
                {
                    return mListOfPdofTimePerType;
                }

                //-------------------------------------------------------------------------------------------------
                // FIXME member function not used
                //void create_unique_adof_list();

                //-------------------------------------------------------------------------------------------------
        };
    }
}

#endif /* SRC_FEM_CL_PDOF_HOST_HPP_ */
