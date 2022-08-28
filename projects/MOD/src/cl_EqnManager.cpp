/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_EqnManager.cpp
 *
 */

#include "cl_EqnManager.hpp" // MOD/src

extern moris::Logger gLogger;

void moris::EqnManager::build_global_system(moris::ArgListEqnMgr_build_global_system & aL)
{
    // Define equation object type
    enum EquationObjectType EqnObjType = STRUC_BAR_1D;

    // Create size of equation object vectors
    moris::Mat< moris::real > ResVec     ( 1,2 );
    moris::Mat< moris::real > JacMat     ( 2,2 );
    moris::Mat< moris::real > MMat       ( 2,2 );
    moris::Mat< moris::real > CMat       ( 2,2 );
    moris::Mat< moris::real > KMat       ( 2,2 );
    moris::Mat< moris::uint > JacMat_col ( 1,4 );
    moris::Mat< moris::uint > JacMat_row ( 1,4 );

    // Fill argument list aL2
    moris::ArgListEqnMgr_get_eqn_obj aL2 = mArgListEqnMgr_get_eqn_obj;
    aL2.EqnObjType = EqnObjType;
    aL2.GlbSolVecU = aL.GlbSolVecU;
    aL2.ResVec     =  ResVec;
    aL2.MMat       =  MMat;
    aL2.CMat       =  CMat;
    aL2.KMat       =  KMat;

    // Loop over equation objects
    for(moris::uint EqnObjId = 0; EqnObjId < mNumberOfEqnObjects; EqnObjId++)
    {
        // Fill current EqnObjId into ArgListEqnMgr_get_eqn_obj
        aL2.EqnObjId   =  EqnObjId;

        // Get equation object size and type
        //myEqnObjPtr.get_eqn_obj_size_and_type(ResVec,JacMat,JacMat_col,JacMat_row,EqnObjType); //PSEUDO CODE

        // Get equations from equation object
        this->get_eqn_obj( aL2);

        // Print outputs for testing purpose
//        aL2.MMat.print("aL2.MMat");
//        aL2.CMat.print("aL2.CMat");
//        aL2.KMat.print("aL2.KMat");
//        aL2.ResVec.print("aL2.ResVec");

        // Build total jacobian through time integrator
        //a2.TimeIntegrator->build_Jacobian (JacMat, JacMat_row, JacMat_col, MMat, CMat, KMat);

#ifdef MORIS_PERFORM_CHECK
        // Perform check of equations
        this->check_eqn_obj(ResVec, JacMat, JacMat_col, JacMat_row,EqnObjId);
#endif

        // Get equation object size and type
        //a2.Solver.sum_into_global_system(ResVec, JacMat,JacMat_col,JacMat_row,ElemId); //PSEUDO CODE
    }
}

void moris::EqnManager::get_eqn_obj( moris::ArgListEqnMgr_get_eqn_obj & aL)
{

    switch (aL.EqnObjType)
    {
    case(STRUC_BAR_1D):
            {

        // Define material properties for structural 1D bar element
        aL.MMat(0,0) =  1;
        aL.MMat(0,1) = -1;
        aL.MMat(1,0) = -1;
        aL.MMat(1,1) =  1;  // For Test purposes
        aL.CMat(0,0) =  1;
        aL.CMat(0,1) = -1;
        aL.CMat(1,0) = -1;
        aL.CMat(1,1) =  1;  // For Test purposes
        aL.KMat(0,0) =  1;
        aL.KMat(0,1) = -1;
        aL.KMat(1,0) = -1;
        aL.KMat(1,1) =  1;  // For Test purposes
        aL.ResVec(0,0) = -1;
        aL.ResVec(0,1) =  1;  // For Test purposes
        break;
            }

    default:
    {
        MORIS_LOG_ERROR ( "Specified equation object type not supported by MORIS");
    }
    }
}

void
moris::EqnManager::check_eqn_obj(moris::Mat<moris::real> &aResVec,
        moris::Mat<moris::real> &aJacMat,
        moris::Mat<moris::uint> &aJacMat_col,
        moris::Mat<moris::uint> &aJacMat_row,
        moris::uint &aEqnObjId)
{

    bool check1 = moris::Math::isfinite( aResVec );
    bool check2 = moris::Math::isfinite( aJacMat );
    bool check3 = moris::Math::isfinite( aJacMat_col );
    bool check4 = moris::Math::isfinite( aJacMat_row );

    if(check1 == false || check2 == false || check3 == false || check4 == false)
    {
        MORIS_LOG_ERROR ( "Equation Object " << aEqnObjId << " has problems");
    }

}

