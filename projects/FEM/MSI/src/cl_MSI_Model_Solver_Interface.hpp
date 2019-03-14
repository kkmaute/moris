/*
 * cl_Model_Solver_Interface.hpp
 *
 *  Created on: Jul 22, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_MANAGER_HPP_
#define SRC_FEM_CL_EQUATION_MANAGER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"

#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Multigrid.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    class Dist_Vector;

    namespace mtk
    {
        class Mesh;
    }
    namespace MSI
    {
        class MSI_Solver_Interface;
        class Model_Solver_Interface
        {
        private:
            moris::Cell< Equation_Object* > & mEquationObjectList;
            Dof_Manager                       mDofMgn;

            mtk::Mesh                       * mMesh;

            Multigrid * mMultigrid = nullptr;

            Param_List< boost::variant< bool, sint, real  > > mMSIParameterList;

            void set_solver_parameters()
            {
                mMSIParameterList.insert( "UX"    ,  1 );
                mMSIParameterList.insert( "UY"    ,  1 );
                mMSIParameterList.insert( "UZ"    ,  1 );
                mMSIParameterList.insert( "TEMP"  ,  1 );
                mMSIParameterList.insert( "L2"    ,  1 );
                mMSIParameterList.insert( "MAPPING_DOF"    ,  1 );
                mMSIParameterList.insert( "LS1",     1 );
                mMSIParameterList.insert( "LS2",     1 );
                mMSIParameterList.insert( "NLSX",    1 );
                mMSIParameterList.insert( "NLSY",    1 );
                mMSIParameterList.insert( "NLSZ",    1 );
                mMSIParameterList.insert( "VX",      1 );
                mMSIParameterList.insert( "VY",      1 );
                mMSIParameterList.insert( "VZ",      1 );
            }


        public:
        Model_Solver_Interface( moris::Cell < Equation_Object* > & aListEqnObj ) : mEquationObjectList( aListEqnObj )
        {
            this->set_solver_parameters();
        };

        /**
         * @brief Model solver interface constructor. This function is tested by the test [MSI_Test][MSI_Test_parallel]
         *
         * @param[in] aListEqnObj   List containing all the equation objects.
         * @param[in] aCommTable    Communication table for adofs.
         *
         */
        Model_Solver_Interface(      moris::Cell < Equation_Object* >                  & aListEqnObj,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs ) : mEquationObjectList( aListEqnObj ),
                                                                                                          mDofMgn( aCommTable, this )
        {
            this->set_solver_parameters();

            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );

            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( aListEqnObj );

            mDofMgn.initialize_pdof_host_list( aListEqnObj );
        };

        Model_Solver_Interface(      moris::Cell < Equation_Object* >                  & aListEqnObj,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs,
                                     mtk::Mesh                                         * aMesh ) : mEquationObjectList( aListEqnObj ),
                                                                                                   mDofMgn( aCommTable, this ),
                                                                                                   mMesh( aMesh )
        {
            this->set_solver_parameters();

            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );

            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( aListEqnObj );

            mDofMgn.initialize_pdof_host_list( aListEqnObj );
        };

//------------------------------------------------------------------------------
        ~Model_Solver_Interface()
        {
            if( mMultigrid != NULL )
            {
                delete mMultigrid;
            }
        };


        void finalize( const bool aUseMultigrid = false )
        {
            mDofMgn.create_adofs();

            mDofMgn.set_pdof_t_matrix();

            for ( Equation_Object* tElement : mEquationObjectList )
            {
                tElement->create_my_pdof_list();

                tElement->create_my_list_of_adof_ids();

                tElement->set_unique_adof_map();

                tElement->set_model_solver_interface_pointer( this );
            }

            if ( aUseMultigrid )
            {
                mMultigrid = new Multigrid( this, mMesh );

                mMultigrid->multigrid_initialize();
            }
        };

//------------------------------------------------------------------------------
        moris::uint get_num_eqn_objs()
        {
            return mEquationObjectList.size();
        };

//------------------------------------------------------------------------------
        Dof_Manager * get_dof_manager()
        {
            return & mDofMgn;
        };

//------------------------------------------------------------------------------

        void
        get_equation_obj_jacobian( const moris::uint      & aEqnObjInd,
                                         Matrix< DDRMat > & aEqnObjMatrix,
                                         Dist_Vector      * aSolutionVector)
        {
            mEquationObjectList( aEqnObjInd )->get_egn_obj_jacobian( aEqnObjMatrix, aSolutionVector );
        };

//------------------------------------------------------------------------------
        void get_equation_obj_residual( const moris::uint      & aEqnObjInd,
                                              Matrix< DDRMat > & aEqnObjRHS,
                                              Dist_Vector      * aSolutionVector )
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_residual( aEqnObjRHS, aSolutionVector  );
        };

//------------------------------------------------------------------------------
        void get_equation_obj_dof_ids( const moris::uint      & aEqnObjInd,
                                             Matrix< DDSMat > & aElementTopology )
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_dof_ids( aElementTopology );
        };

//------------------------------------------------------------------------------
        void read_multigrid_maps( const moris::uint               aLevel,
                                  const moris::Matrix< DDSMat > & aExtFineIndices,
                                  const moris::sint               aTypeTimeIdentifier,
                                        moris::Matrix< DDSMat > & aInternalFineIndices)
        {
            mMultigrid->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
        };

//------------------------------------------------------------------------------
        const moris::Cell< Matrix< DDUMat > > & get_lists_of_ext_index_multigrid()
        {
            return mMultigrid->get_lists_of_ext_index_multigrid();
        };

        const moris::Cell< Matrix< DDSMat > > & get_lists_of_multigrid_identifiers()
        {
            return mMultigrid->get_lists_of_multigrid_identifiers();
        };

        const moris::Cell< moris::Cell< Matrix< DDSMat > > > & get_multigrid_map( )
        {
            return mMultigrid->get_multigrid_map();
        };

        const Matrix< DDUMat > & get_number_remaining_dofs()
        {
            return mMultigrid->get_number_remaining_dofs();
        };

//------------------------------------------------------------------------------
        mtk::Mesh * get_mesh_pointer_for_multigrid( )
        {
            return mMesh;
        };

        Multigrid * get_MSI_multigrid_Pointer( )
        {
            return mMultigrid;
        };

        boost::variant< bool, sint, real > & set_param( const std::string & aKey )
        {
            return mMSIParameterList( aKey );
        }

        moris::sint get_adof_order_for_type( moris::uint aDofType )
        {
           // Get dof type enum
           enum Dof_Type tDofType = mDofMgn.get_dof_type_enum( aDofType );

           if      ( tDofType == Dof_Type::UX )          { return mMSIParameterList.get< moris::sint >( "UX" ); }
           else if ( tDofType == Dof_Type::UY )          { return mMSIParameterList.get< moris::sint >( "UY" ); }
           else if ( tDofType == Dof_Type::UZ )          { return mMSIParameterList.get< moris::sint >( "UZ" ); }
           else if ( tDofType == Dof_Type::TEMP )        { return mMSIParameterList.get< moris::sint >( "TEMP" ); }
           else if ( tDofType == Dof_Type::L2 )          { return mMSIParameterList.get< moris::sint >( "L2" ); }
           else if ( tDofType == Dof_Type::MAPPING_DOF ) { return mMSIParameterList.get< moris::sint >( "MAPPING_DOF" ); }
           else if ( tDofType == Dof_Type::LS1 )         { return mMSIParameterList.get< moris::sint >( "LS1" ); }
           else if ( tDofType == Dof_Type::LS2 )         { return mMSIParameterList.get< moris::sint >( "LS2" ); }
           else if ( tDofType == Dof_Type::NLSX )        { return mMSIParameterList.get< moris::sint >( "NLSX" ); }
           else if ( tDofType == Dof_Type::NLSY )        { return mMSIParameterList.get< moris::sint >( "NLSY" ); }
           else if ( tDofType == Dof_Type::NLSZ )        { return mMSIParameterList.get< moris::sint >( "NLSZ" ); }
           else if ( tDofType == Dof_Type::VX )          { return mMSIParameterList.get< moris::sint >( "VX" ); }
           else if ( tDofType == Dof_Type::VY )          { return mMSIParameterList.get< moris::sint >( "VY" ); }
           else if ( tDofType == Dof_Type::VZ )          { return mMSIParameterList.get< moris::sint >( "VZ" ); }

           else
           {
               MORIS_ERROR( false, "Model_Solver_Interface::get_adof_order_for_type(): Dof type does not exist. Check dof type enums");
               return 0;
           }
        }


    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_MANAGER_HPP_ */
