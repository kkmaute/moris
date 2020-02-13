/*
 * cl_MSI_Model_Solver_Interface.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Model_Solver_Interface.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_MSI_Equation_Model.hpp"

namespace moris
{
    namespace MSI
    {

//------------------------------------------------------------------------------

    Model_Solver_Interface::Model_Solver_Interface(      ParameterList                                       aMSIParameterList,
                                                         std::shared_ptr< MSI::Equation_Model >            & aEquationModel,
                                                   const Matrix< IdMat >                                   & aCommTable,
                                                   const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                                                   const moris::uint                                         aNumMaxAdofs,
                                                         mtk::Mesh                                         * aMesh ) : mMSIParameterList( aMSIParameterList ),
                                                                                                                       mEquationBlocks( aEquationModel->get_equation_sets() ),
                                                                                                                       mDofMgn( aCommTable, this ),
                                                                                                                       mMesh( aMesh )
    {
        this->create_equation_object_list();

        mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );

        mDofMgn.set_max_num_adofs( aNumMaxAdofs );

        mDofMgn.initialize_pdof_type_list( mEquationBlocks );
    }

//------------------------------------------------------------------------------

    moris::sint Model_Solver_Interface::get_adof_index_for_type( moris::uint aDofType )
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
       else if ( tDofType == Dof_Type::P )           { return mMSIParameterList.get< moris::sint >( "P" ); }

       else
       {
           MORIS_ERROR( false, "Model_Solver_Interface::get_adof_index_for_type(): Dof type does not exist. Check dof type enums");
           return 0;
       }
    }

    //------------------------------------------------------------------------------

    moris::sint Model_Solver_Interface::get_max_adof_index()
    {
        // max BSpline mesh index. However, not limited to BSpline meshes.
        sint tMaxAdofIndex = 0;

        uint tNumDofTypes = mDofMgn.get_num_dof_types();

        // loop over all used Ddf type indices
        for( uint Ik = 0; Ik< tNumDofTypes; Ik++)
        {
            sint tIndex = this->get_adof_index_for_type( Ik );

            tMaxAdofIndex = std::max( tIndex, tMaxAdofIndex );
        }

        return tMaxAdofIndex+1;
    }

    }
}
