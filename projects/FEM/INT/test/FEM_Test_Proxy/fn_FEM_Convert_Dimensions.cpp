/*
 * fn_FEM_Convert_Dimensions.cpp
 *
 *  Created on: Jun 23, 2021
 *      Author: wunsch
 */

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        /* The below functions assume the 2D element to be a QUAD9 - TIME2 element 
        * with 18 Nodes and the 1D element to be a LINE3 - TIME2 element with 6 DoFs
        * 
        * Ordering Scheme:
        * 
        * Timelevels:      t_1            t_2
        * 
        *        1D:   1 -- 3 -- 2    4 --  6 -- 5
        * 
        *              4 -- 7 -- 3   13 -- 16 -- 12
        *              |         |    |          |
        *        2D:   8    9    6   17    18    15
        *              |         |    |          |
        *              1 -- 5 -- 2   10 -- 14 -- 11
        */

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_residual_2D_to_1D_quadratic(
                const Matrix< DDRMat > & aRes2D )
        {
            // check size of the residual inputed
            MORIS_ERROR( ( aRes2D.n_rows() == 18 ) and ( aRes2D.n_cols() == 1 ),
                    "convert_residual_2D_to_1D_quadratic() : 2D residual of wrong size, must be 18 x 1 " );

            // convert the 2D residual to 1D
            Matrix< DDRMat > tRes1D = { 
                    { aRes2D(  0 ) + aRes2D(  3 ) + aRes2D(  7 ) },
                    { aRes2D(  1 ) + aRes2D(  2 ) + aRes2D(  5 ) },
                    { aRes2D(  4 ) + aRes2D(  6 ) + aRes2D(  8 ) },
                    { aRes2D(  9 ) + aRes2D( 12 ) + aRes2D( 16 ) },
                    { aRes2D( 10 ) + aRes2D( 11 ) + aRes2D( 14 ) },
                    { aRes2D( 13 ) + aRes2D( 15 ) + aRes2D( 17 ) } 
                    };

            // return converted residual vector
            return tRes1D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_multi_residual_2D_to_1D_quadratic(
                const Matrix< DDRMat > & aMultiRes2D,
                const uint               aNumStateVars )
        {
            // check size of the residual inputed
            MORIS_ERROR( ( aMultiRes2D.n_rows() == 18 * aNumStateVars ) and ( aMultiRes2D.n_cols() == 1 ) and ( aNumStateVars > 0 ),
                    "convert_multi_residual_2D_to_1D_quadratic() : 2D residual of wrong size, must be (#StateVars * 18) x 1 " );

            // initialize converted residual vector
            Matrix< DDRMat > tMultiRes1D( aNumStateVars * 6, 1, 0.0 );

            // convert 2D res vec for each state variable and place into 1D residual vector
            for ( uint iVar = 0; iVar < aNumStateVars; iVar++ )
            {
                // get section of state variable vector corresponding to current state var
                Matrix< DDRMat > tSingleVarRes2D = aMultiRes2D( { iVar * 18, ( iVar + 1 ) * 18 - 1 } );

                // convert 2D residual to 1D
                Matrix< DDRMat > tSingleVarRes1D = convert_residual_2D_to_1D_quadratic( tSingleVarRes2D );

                // fill converted 1D residual vector
                tMultiRes1D( { iVar * 6, ( iVar + 1 ) * 6 - 1 } ) = tSingleVarRes1D.matrix_data();
            }

            // return converted residual vector
            return tMultiRes1D;
        }

        //------------------------------------------------------------------------------  

        Matrix< DDRMat > convert_comp_flow_residual_2D_to_1D_quadratic( 
                const Matrix< DDRMat > & aCFRes2D )
        {
            // check size of the residual inputed
            MORIS_ERROR( ( aCFRes2D.n_rows() == 18 * 4 ) and ( aCFRes2D.n_cols() == 1 ),
                    "convert_comp_flow_residual_2D_to_1D_quadratic() : 2D residual of wrong size, must be 72 x 1 " );

            // convert multi-state-var residual to 1D
            Matrix< DDRMat > tCFRes1D = convert_multi_residual_2D_to_1D_quadratic( aCFRes2D, 4 );

            // cut out UY part
            tCFRes1D( { 12, 17 } ) = tCFRes1D( { 18, 23 } );
            tCFRes1D = tCFRes1D( { 0, 17 } );

            // return 
            return tCFRes1D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_Jacobian_2D_to_1D_quadratic (
                const Matrix< DDRMat > & aJac2D )
        {
            // check size of the Jacobian inputed
            MORIS_ERROR( ( aJac2D.n_rows() == 18 ) and ( aJac2D.n_cols() == 18 ),
                    "convert_Jacobian_2D_to_1D_quadratic() : 2D jacobian of wrong size, must be 18 x 18 " );

            // // lumping matrix defining mapping the entries of the 1D residual / jacobian go to the 2D entries
            // Matrix< DDRMat > tLumpMat = {
            //         {  0,  3,  7 },
            //         {  1,  2,  5 },
            //         {  4,  6,  8 },
            //         {  9, 12, 16 },
            //         { 10, 11, 14 },
            //         { 13, 15, 17 } };

            // inverse lumping vector mapping the entries of the 2D residual / jacobian go to the 1D entries
            Matrix< DDRMat > tInvLumpMat = { { 0, 1, 1, 0, 2, 1, 2, 0, 2, 3, 4, 4, 3, 5, 4, 5, 3, 5 } };

            // initialize 1D jacobian
            Matrix< DDRMat > tJac1D( 6, 6, 0.0 );

            // loop over the rows of the 2D jacobian, lumping together the rows of the 1D matrix
            for ( uint iRow2D = 0; iRow2D < tInvLumpMat.length(); iRow2D++ )
            {
                // get the corresponding 1D row
                uint tRow1D = tInvLumpMat( iRow2D );

                // loop over the rows of the 2D jacobian, lumping together the rows of the 1D matrix
                for ( uint iCol2D = 0; iCol2D < tInvLumpMat.length(); iCol2D++ )
                {
                    // get the corresponding 1D column
                    uint tCol1D = tInvLumpMat( iCol2D );

                    // add entries together
                    tJac1D( tRow1D, tCol1D ) = tJac1D( tRow1D, tCol1D ) + aJac2D( iRow2D, iCol2D );
                }
            }
            
            // return 1D jacobian
            return tJac1D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_multi_Jacobian_2D_to_1D_quadratic (
                const Matrix< DDRMat > & aJac2D,
                const uint               aNumStateVars )
        {
            // check size of the Jacobian inputed
            MORIS_ERROR( ( aJac2D.n_rows() == 18 * aNumStateVars ) and ( aJac2D.n_cols() == 18 * aNumStateVars ) and ( aNumStateVars > 0 ),
                    "convert_multi_Jacobian_2D_to_1D_quadratic() : 2D Jacobian of wrong size, must be (#StateVars * 18) x (#StateVars * 18) " );

            // initialize multi-state-var 1D jacobian
            Matrix< DDRMat > tMultiJac1D( aNumStateVars * 6, aNumStateVars * 6, 0.0 );

            // loop over each state variable block in the jacobian
            for ( uint iRowVarBlock = 0; iRowVarBlock < aNumStateVars; iRowVarBlock++ )
            {
                for ( uint iColVarBlock = 0; iColVarBlock < aNumStateVars; iColVarBlock++ )
                {
                    // extract block from 2D jacobian
                    Matrix< DDRMat > tJac2DBlock = aJac2D( 
                            { iRowVarBlock * 18, ( iRowVarBlock + 1 ) * 18 - 1 }, 
                            { iColVarBlock * 18, ( iColVarBlock + 1 ) * 18 - 1 } );

                    // convert block to 1D
                    Matrix< DDRMat > tJac1DBlock = convert_Jacobian_2D_to_1D_quadratic( tJac2DBlock );

                    // put 1D block onto 1D Jacobian
                    tMultiJac1D( 
                            { iRowVarBlock * 6, ( iRowVarBlock + 1 ) * 6 - 1 }, 
                            { iColVarBlock * 6, ( iColVarBlock + 1 ) * 6 - 1 } ) = tJac1DBlock.matrix_data();
                }
            }

            // return 1D jacobian
            return tMultiJac1D;
        }

        //------------------------------------------------------------------------------  

        Matrix< DDRMat > convert_comp_flow_jacobian_2D_to_1D_quadratic( 
                const Matrix< DDRMat > & aCFJac2D )
        {
            // check size of the Jacobian inputed
            MORIS_ERROR( ( aCFJac2D.n_rows() == 18 * 4 ) and ( aCFJac2D.n_cols() == 18 * 4 ),
                    "convert_comp_flow_Jacobian_2D_to_1D_quadratic() : 2D Jacobian of wrong size, must be 72 x 72 " );

            // convert multi-state-var residual to 1D
            Matrix< DDRMat > tCFJac1D = convert_multi_Jacobian_2D_to_1D_quadratic( aCFJac2D, 4 );

            // cut out UY part
            tCFJac1D( {  0, 11 }, { 12, 17 } ) = tCFJac1D( {  0, 11 }, { 18, 23 } );
            tCFJac1D( { 12, 17 }, {  0, 11 } ) = tCFJac1D( { 18, 23 }, {  0, 11 } );
            tCFJac1D( { 12, 17 }, { 12, 17 } ) = tCFJac1D( { 18, 23 }, { 18, 23 } );
            tCFJac1D = tCFJac1D( { 0, 17 }, { 0, 17 } );

            // return 
            return tCFJac1D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_DoF_vector_1D_to_2D_quadratic (
                const Matrix< DDRMat > & aDofVec1D )
        {
            // check size of the DoF vector inputed
            MORIS_ERROR( aDofVec1D.n_rows() == 6 && aDofVec1D.n_cols() == 1,
                    "convert_DoF_vector_1D_to_2D_quadratic() : 1D DoF vector of wrong size, must be 6 x 1 " );

            // DoFs into 2D DoF vector
            Matrix< DDRMat > tDofVec2D = {
                    { aDofVec1D( 0 ) },
                    { aDofVec1D( 1 ) },
                    { aDofVec1D( 1 ) },
                    { aDofVec1D( 0 ) },
                    { aDofVec1D( 2 ) },
                    { aDofVec1D( 1 ) },
                    { aDofVec1D( 2 ) },
                    { aDofVec1D( 0 ) },
                    { aDofVec1D( 2 ) },
                    { aDofVec1D( 0 + 3 ) },
                    { aDofVec1D( 1 + 3 ) },
                    { aDofVec1D( 1 + 3 ) },
                    { aDofVec1D( 0 + 3 ) },
                    { aDofVec1D( 2 + 3 ) },
                    { aDofVec1D( 1 + 3 ) },
                    { aDofVec1D( 2 + 3 ) },
                    { aDofVec1D( 0 + 3 ) },
                    { aDofVec1D( 2 + 3 ) } 
                    };

            // return converted DoF vector
            return tDofVec2D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_multi_DoF_vector_1D_to_2D_quadratic (
                const Matrix< DDRMat > & aDofVec1D,
                const uint               aNumStateVars )
        {
            // check size of the DoF vector inputed
            MORIS_ERROR( ( aDofVec1D.n_rows() == 6 * aNumStateVars ) and ( aDofVec1D.n_cols() == 1 ),
                    "convert_multi_DoF_vector_1D_to_2D_quadratic() : 1D DoF vector of wrong size, must be (#StateVars * 6) x 1 " );

            // initialize converted dof vector
            Matrix< DDRMat > tMultiDofVec2D( aNumStateVars * 18, 1, 0.0 );

            // convert 2D dof vec for each state variable and place into 2D dof vector
            for ( uint iVar = 0; iVar < aNumStateVars; iVar++ )
            {
                // get section of 1D dof vector corresponding to current state var
                Matrix< DDRMat > tSingleVarDofVec1D = aDofVec1D( { iVar * 6, ( iVar + 1 ) * 6 - 1 } );

                // convert 1D dof vec to 2D
                Matrix< DDRMat > tSingleVarDofVec2D = convert_DoF_vector_1D_to_2D_quadratic( tSingleVarDofVec1D );

                // fill converted 1D residual vector
                tMultiDofVec2D( { iVar * 18, ( iVar + 1 ) * 18 - 1 } ) = tSingleVarDofVec2D.matrix_data();
            }

            // return converted DoF vector
            return tMultiDofVec2D;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > convert_comp_flow_DoF_vector_1D_to_2D_quadratic (
                const Matrix< DDRMat > & aCFDofVec1D )
        {
            // check size of the DoF vector inputed
            MORIS_ERROR( ( aCFDofVec1D.n_rows() == 6 * 3 ) and ( aCFDofVec1D.n_cols() == 1 ),
                    "convert_comp_flow_DoF_vector_1D_to_2D_quadratic() : 1D DoF vector of wrong size, must be 18 x 1 " );

            // initialize converted dof vector
            Matrix< DDRMat > tCFDofVec2D( 72, 1, 0.0 );

            // convert dof vector
            Matrix< DDRMat > tCFDofVec2Dtemp = convert_multi_DoF_vector_1D_to_2D_quadratic( aCFDofVec1D, 3 );

            // add zeros for UY to converted Dof vector
            tCFDofVec2D( {  0, 35 } ) = tCFDofVec2Dtemp( {  0, 35 } );
            tCFDofVec2D( { 54, 71 } ) = tCFDofVec2Dtemp( { 36, 53 } );

            // return converted DoF vector
            return tCFDofVec2D;
        }

        //------------------------------------------------------------------------------

    } // end namespace fem
} // end namespace moris