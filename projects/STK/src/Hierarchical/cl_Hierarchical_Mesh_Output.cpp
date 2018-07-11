/*
 * cl_Hierarchical_Mesh_Output.cpp
 *
 *  Created on: Jan 8, 2018
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Output.hpp" // STK/src/Hierarchical
using namespace moris;

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_element_dof_connectivity_ascii(
        const uint & aDim,
        const uint & aPolynomial,
        const Mat<uint> & aNumElements,
        const Mat<uint> & aFetopo,
        const Mat<uint> & aElementListOnProc,
        const BoostBitset & aBasisActive,
        const bool & aRefinement,
        const bool & aTruncatedBsplines,
        Cell<Mat<real>> & aTMatrixParentChild,
        Mat<uint> & aBsplineBasisList,
        const map< uint, uint > & aBasisListMap )
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();

    // Write Dof connectivity for FEMDOC
    Mat<uint> tOrder(pow((aPolynomial)+1,(aDim)),1);
    uint tVar = 0;
    Mat<real> tHelp;
    Mat<uint> tHelpList;
    std::ofstream DofNodeConnecitivity;
    if(  aRefinement == 1)
    {
        if( tProcSize == 1)
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data");
        }
        else
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data.org");
        }
        else
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data.org_" + std::to_string(tProcRank));
        }
    }

    // FIXME: move this into create_maps
    //Create the map for the consecutive numbering
    for(uint i=0; i < aBsplineBasisList.length(); i++)
    {
        aBsplineBasisList(i) = i+1;
    }

    Mat<uint> tIdField;
    Mat<real> tTMatrix;
    for(uint i=0; i<aElementListOnProc.length(); i++)
    {
        if( aTruncatedBsplines == false)
        {
            mTMatrix.give_Tmatrix_and_IdField_Reorder(
                    aElementListOnProc(i),
                    aDim,
                    aPolynomial,
                    aNumElements,
                    aBasisActive,
                    tTMatrix,
                    tIdField,
                    aTMatrixParentChild);
        }
        else
        {
            mTMatrix.give_Truncated_Tmatrix_and_IdField_Reorder(
                    aElementListOnProc(i),
                    aDim,
                    aPolynomial,
                    aNumElements,
                    aBasisActive,
                    tTMatrix,
                    tIdField,
                    aTMatrixParentChild);
        }
        tHelpList.set_size(tTMatrix.n_rows(),1,0);
        tVar = 0;
        for(uint j=0; j<tIdField.length(); j++)
        {
            tHelp = tTMatrix.row(j);
            if( sum( tHelp ) != 0.0 )
            {
                tHelpList(tVar) = j;
                tVar++;
            }
        }
        tHelpList.resize(tVar,1);
        DofNodeConnecitivity << aElementListOnProc(i) << "\n";
        DofNodeConnecitivity <<  tHelpList.length() << " " << pow((aPolynomial+1),aDim) << "\n";

        MORIS_ASSERT( tIdField.length() >= pow((aPolynomial+1),aDim), "There are less then the minimum number of DOFS active ." );

        for(uint j=0; j<tHelpList.length(); j++)
        {
            for(uint k=0; k<pow((aPolynomial+1),aDim); k++)
            {
                DofNodeConnecitivity << tTMatrix(tHelpList(j),k) << " ";
            }
            DofNodeConnecitivity << "\n";
        }
        for(uint j=0; j<tHelpList.length(); j++)
        {
            DofNodeConnecitivity << aBasisListMap.find(tIdField(tHelpList(j)))  << "\n";
        }
    }
    DofNodeConnecitivity.close();
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_element_dof_connectivity_binary(
        const uint              & aDim,
        const uint              & aPolynomial,
        const Mat< uint >       & aNumElements,
        const Mat< uint >       & aFetopo,
        const Mat< uint >       & aElementListOnProc,
        const BoostBitset       & aBasisActive,
        const bool              & aRefinement,
        const bool              & aTruncatedBsplines,
        Cell< Mat< real > >     & aTMatrixParentChild,
        Mat<uint>               & aBsplineBasisList,
        const map< uint, uint > & aBasisListMap )
{

    // FIXME: move this into create_maps
    //Create the map for the consecutive numbering
    for(uint i=0; i < aBsplineBasisList.length(); i++)
    {
        aBsplineBasisList(i) = i+1;
    }

    // get number of elements
    uint tNumberOfElements = aElementListOnProc.length();

    // T-Matrices and indices for current element
    Cell < Mat< real > > tTMatrices;
    Cell < Mat< uint > > tIdFields;
    Cell < Mat< uint > > tIndices;

    // temporary variable containing the row of a T-Matrix
    Mat<real> tRow;

    // reserve memory
    tTMatrices.resize( tNumberOfElements );
    tIdFields.resize( tNumberOfElements );
    tIndices.resize( tNumberOfElements );

    // number of Lagrange nodes per element
    uint tNumberOfLagrangeNodesPerElement = pow( ( aPolynomial+1 ), aDim );

    // loop over all elements and calculate T-Matrix and ID field
    if( aTruncatedBsplines == false)
    {
        // loop over all elements
        for( uint e=0; e<tNumberOfElements; ++e )
        {
            Mat < uint > tIdField;
            Mat < real > tTMatrix;

            // calculate T-Matrix and save into cell
            mTMatrix.give_Tmatrix_and_IdField_Reorder(
                    aElementListOnProc( e ),
                    aDim,
                    aPolynomial,
                    aNumElements,
                    aBasisActive,
                    tTMatrix,
                    tIdFields( e ),
                    aTMatrixParentChild);

            // check length of ID field
            MORIS_ASSERT( tIdFields( e ).length() >= tNumberOfLagrangeNodesPerElement ,
                    "There are less then the minimum number of DOFs active ." );
            // create help vector
            tIndices(e).set_size(tTMatrices(e).n_rows(),1,0);
            uint tCount = 0;
            for (uint j=0; j<tIdFields( e ).length(); ++j )
            {
                tRow = tTMatrices( e ).row( j );
                if( sum( tRow ) != 0.0 )
                {
                    tIndices( e )( tCount++ ) = j;
                }
            }

            // resize help index
            tIndices(e).resize( tCount, 1 );
        }
    }
    else
    {
        // loop over all elements
        for( uint e=0; e<tNumberOfElements; ++e )
        {
            // calculate T-Matrix and save into cell
            mTMatrix.give_Truncated_Tmatrix_and_IdField_Reorder(
                    aElementListOnProc( e ),
                    aDim,
                    aPolynomial,
                    aNumElements,
                    aBasisActive,
                    tTMatrices( e ),
                    tIdFields( e ),
                    aTMatrixParentChild);

            // check length of ID field
            MORIS_ASSERT( tIdFields( e ).length() >= tNumberOfLagrangeNodesPerElement ,
                    "There are less then the minimum number of DOFs active ." );

            // create help vector
            tIndices(e).set_size( tTMatrices( e ).n_rows(),1,0);
            uint tCount = 0;
            for (uint j=0; j<tIdFields( e ).length(); ++j )
            {
                tRow = tTMatrices(e).row(j);
                if( sum( tRow ) != 0.0 )
                {
                    tIndices( e )( tCount++ ) = j;
                }
            }

            // resize help index
            tIndices(e).resize( tCount, 1 );
        }
    }
    // count samples needed for output matrix
    uint tCount = 0;

    // loop over all elements
    for( uint e=0; e<tNumberOfElements; ++e )
    {
        // increment one for element ID
        ++tCount;

        // increment two for T-Matrix size
        tCount += 2;

        // increment about size of matrix
        tCount += tIndices( e ).length() * tNumberOfLagrangeNodesPerElement;

        // increment about size of indices
        tCount += tIndices( e ).length();

    }

    // assign memory of output matrix
    Mat < real > tDofConnectivity( tCount, 1 );

    // reset counter
    tCount = 0;

    // loop over all elements
    for( uint e=0; e<tNumberOfElements; ++e )
    {
        // save element ID
        tDofConnectivity( tCount++ ) = ( real) aElementListOnProc( e );

        // size of T-Matrix
        tDofConnectivity( tCount++ ) = ( real ) tIndices( e ).length();

        // number of Lagrange DOFs
        tDofConnectivity( tCount++ ) = ( real ) tNumberOfLagrangeNodesPerElement ;

        // save T-Matrix in row major order
        for ( uint j=0; j<tIndices( e ).length(); ++j)
        {
            for ( uint i=0; i< tNumberOfLagrangeNodesPerElement; ++i)
            {
                tDofConnectivity( tCount++ ) = tTMatrices( e )( tIndices( e )( j ), i );
            }
        }

        // save B-Spline DOFs
        for(uint j=0; j<tIndices( e ).length(); j++)
        {

            tDofConnectivity( tCount++ )
                                = ( real ) aBasisListMap.find( tIdFields( e )( tIndices( e )( j ) ) );

        }

    }

    MORIS_ASSERT ( tDofConnectivity.length() == tCount,
                   "Size of tDofConnectivity is wrong" );

    // get name for output file
    std::string tOutputFile;

    // pick name
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tOutputFile = "ElementalDofNodeConnecitivity.data";
        }
        else
        {
            tOutputFile = "ElementalDofNodeConnecitivity.data_" + std::to_string( par_rank() );
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tOutputFile = "ElementalDofNodeConnecitivity.data.org";
        }
        else
        {
            tOutputFile = "ElementalDofNodeConnecitivity.data.org_" + std::to_string( par_rank() );
        }
    }

    // write output vector to file
    save_vector_to_binary_file( tOutputFile, tDofConnectivity );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_node_dof_connectivity_ascii(
        const Mat<uint>         & aLagrangeToBSplineMap,
        const Mat<uint>         & aIdFieldField,
        const Mat<real>         & aTMatrixField,
        const map< uint, uint > & aBasisListMap,
        const bool              & aRefinement )
{
    // Write Dof connectivity for FEMDOC
    std::ofstream DOFVariables;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            DOFVariables.open ("NodalDofNodeConnecitivity.data");
        }
        else
        {
            DOFVariables.open ("NodalDofNodeConnecitivity.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            DOFVariables.open ("NodalDofNodeConnecitivity.data.org");
        }
        else
        {
            DOFVariables.open ("NodalDofNodeConnecitivity.data.org_" + std::to_string(par_rank()));
        }
    }

    // loop over all active lagrange basis on proc
    for(uint i=0; i<aLagrangeToBSplineMap.length(); i++)
    {
        // write global lagrange basis ID and 2 * number of associated B-Spline basis
        DOFVariables << aLagrangeToBSplineMap(i) << " " << 2*aTMatrixField(i,0) << "\n";

        for(uint j = 0; j< aIdFieldField(i,0); j++)
        {
            DOFVariables << aBasisListMap.find(aIdFieldField(i,j+1)) << " ";
        }

        for(uint j = 0; j< aIdFieldField(i,0); j++)
        {
            DOFVariables << aTMatrixField(i,j+1) << " ";
        }
        DOFVariables << "\n";
    }
    DOFVariables.close();
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_node_dof_connectivity_binary(
        const Mat<uint>         & aLagrangeToBSplineMap,
        const Mat<uint>         & aIdFieldField,
        const Mat<real>         & aTMatrixField,
        const map< uint, uint > & aBasisListMap,
        const bool              & aRefinement)
{
    // calculate number of samples
    uint tCount = 0;

    for( uint i=0; i<aLagrangeToBSplineMap.length(); ++i )
    {
        // increment is
        // 1 sample for node id
        // 1 sample for n samples for node
        // n samples
        tCount += 2 + 2*aTMatrixField( i, 0 );
    }

    // assign memory of output matrix
    Mat < real > tDofConnectivity( tCount, 1 );

    // reset counter
    tCount = 0;

    // loop over all active lagrange basis on proc
    for(uint i=0; i<aLagrangeToBSplineMap.length(); i++)
    {
        // write global Lagrange basis ID and 2 * number of associated B-Spline basis
        tDofConnectivity( tCount++ ) = aLagrangeToBSplineMap( i );

        // write 2 * number of associated B-Spline basis
        tDofConnectivity( tCount++ ) = 2*aTMatrixField(i,0);

        // write associated B-Spline basis in new numbering
        for(uint j = 0; j< aIdFieldField(i,0); j++)
        {
            tDofConnectivity( tCount++ ) = aBasisListMap.find(aIdFieldField(i,j+1));
        }

        // write entries of T-Matrix
        for(uint j = 0; j< aIdFieldField(i,0); j++)
        {
            tDofConnectivity( tCount++ ) = aTMatrixField(i,j+1);
        }
    }

    // make sure that everything is correct
    MORIS_ASSERT ( tDofConnectivity.length() == tCount,
                       "Size of tDofConnectivity is wrong" );

    // get name for output file
    std::string tOutputFile;

    // pick name for output file
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tOutputFile = "NodalDofNodeConnecitivity.data";
        }
        else
        {
            tOutputFile = "NodalDofNodeConnecitivity.data_" + std::to_string(par_rank());
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tOutputFile = "NodalDofNodeConnecitivity.data.org";
        }
        else
        {
            tOutputFile = "NodalDofNodeConnecitivity.data.org_" + std::to_string(par_rank());
        }
    }

    // write output vector to file
    save_vector_to_binary_file( tOutputFile, tDofConnectivity );

}

// -----------------------------------------------------------------------------
void Hierarchical_Mesh_Output::output_design_variables_for_moris_ascii(
        const Mat<uint> & aLagrangeToBSplineMap,
        const Mat<uint> & aIdFieldFieldDesign,
        const Mat<real> & aTMatrixFieldDesign,
        const map<uint, uint> & aDesignBSplineListMap,
        const bool & aRefinement,
        const bool & aSymmetry
)

{
    // Write Dof connectivity for MORIS

    std::ofstream DesignVariablesMoris;
    // vector is almost identical except for two sizes in the beginning
    if(  aRefinement > 0)
    {
        if( par_size() == 1)
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris_new.data");
        }
        else
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris.data.org");
        }
        else
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris.data.org_" + std::to_string(par_rank()));
        }
    }
    DesignVariablesMoris << aIdFieldFieldDesign.n_rows() << "\n";
    DesignVariablesMoris << aIdFieldFieldDesign.n_cols() << "\n";
    if ( ! aSymmetry  )
    {
        // standard case

        // loop over all lagrange nodes
        for(uint i=0; i<aLagrangeToBSplineMap.length(); i++)
        {

            DesignVariablesMoris << aLagrangeToBSplineMap(i) << "\n";
            DesignVariablesMoris << 2*aTMatrixFieldDesign(i,0) << "\n";
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                DesignVariablesMoris << aDesignBSplineListMap.find(aIdFieldFieldDesign(i,j+1)) << "\n";
            }
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                DesignVariablesMoris << aTMatrixFieldDesign(i,j+1) << "\n";
            }
        }
    }
    else
    {
        // symmetric case
        // loop over all Lagrange nodes
        for(uint i=0; i<aLagrangeToBSplineMap.length(); i++ )
        {
            // calculate number of samples
            uint tNumberOfSamples = aTMatrixFieldDesign( i, 0 );

            // allocate memory
            Mat< uint > tIdField( tNumberOfSamples , 1 );
            Mat< real > tTMatrix ( tNumberOfSamples, 1 );

            // copy values from input
            for( uint j=0; j<tNumberOfSamples; ++j )
            {
                tIdField( j ) =  aDesignBSplineListMap.find( aIdFieldFieldDesign(i,j+1) );
                tTMatrix( j ) =  aTMatrixFieldDesign(i,j+1);
            }

            // make unique
            this-> make_design_variables_unique( tIdField, tTMatrix );

            // readjust number of samples
            tNumberOfSamples = tIdField.length();

            // write output
            DesignVariablesMoris << aLagrangeToBSplineMap(i) << "\n";
            DesignVariablesMoris << 2*tNumberOfSamples << "\n";
            for(uint j = 0; j< tNumberOfSamples; ++j )
            {
                DesignVariablesMoris << tIdField( j ) << "\n";
            }
            for(uint j = 0; j< tNumberOfSamples; ++j )
            {
                DesignVariablesMoris << tTMatrix( j ) << "\n";
            }
        }
    }

    DesignVariablesMoris.close();
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_design_variables_for_femdoc_ascii(
        const Mat<uint> & aLagrangeToBSplineMap,
        const Mat<uint> & aIdFieldFieldDesign,
        const Mat<real> & aTMatrixFieldDesign,
        const map<uint, uint> & aDesignBSplineListMap,
        const bool & aRefinement,
        const bool & aSymmetry
        )

{

    // Write Dof connectivity for FEMDOC
    std::ofstream DesignVariables;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            DesignVariables.open ("DesignVariables_new.data");
        }
        else
        {
            DesignVariables.open ("DesignVariables_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            DesignVariables.open ("DesignVariables.data.org");
        }
        else
        {
            DesignVariables.open ("DesignVariables.data.org_" + std::to_string(par_rank()));
        }
    }

    if ( ! aSymmetry  )
    {
        // standard case
        for(uint i=0; i < aLagrangeToBSplineMap.length(); i++)
        {
            // DesignVariables << aNodalLocaltoGlobal(i) << " " << 2*aTMatrixFieldDesign(i,0) << "\n";

            DesignVariables << aLagrangeToBSplineMap(i) << " " << 2*aTMatrixFieldDesign(i,0) << "\n";

            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                DesignVariables << aDesignBSplineListMap.find(aIdFieldFieldDesign(i,j+1)) << " ";
            }
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                DesignVariables << aTMatrixFieldDesign(i,j+1) << " ";
            }
            DesignVariables << "\n";
        }
    }
    else
    {
        // symmetric case

        // loop over all indices
        for(uint i=0; i < aLagrangeToBSplineMap.length(); i++)
        {
            // calculate number of samples
            uint tNumberOfSamples = aTMatrixFieldDesign(i,0);

            // allocate memory
            Mat< uint > tIdField( tNumberOfSamples , 1 );
            Mat< real > tTMatrix ( tNumberOfSamples, 1 );

            // copy values from input
            for( uint j=0; j<tNumberOfSamples; ++j )
            {
                tIdField( j ) =  aDesignBSplineListMap.find( aIdFieldFieldDesign(i,j+1) );
                tTMatrix( j ) =  aTMatrixFieldDesign(i,j+1);
            }

            // combine symmetric entries
            this-> make_design_variables_unique( tIdField, tTMatrix );

            // readjust number of samples
            tNumberOfSamples = tIdField.length();

            // write to output
            DesignVariables << aLagrangeToBSplineMap(i) << " " << 2*tNumberOfSamples << "\n";
            for(uint j = 0; j< tNumberOfSamples; ++j )
            {
                DesignVariables << tIdField( j ) << " ";
            }
            for(uint j = 0; j< tNumberOfSamples; ++j)
            {
                DesignVariables << tTMatrix( j ) << " ";
            }

            // finish line
            DesignVariables << "\n";
        }
    }
    DesignVariables.close();
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_design_variables_for_femdoc_binary(
        const Mat<uint>       & aLagrangeToBSplineMap,
        const Mat<uint>       & aIdFieldFieldDesign,
        const Mat<real>       & aTMatrixFieldDesign,
        const map<uint, uint> & aDesignBSplineListMap,
        const bool            & aRefinement,
        const bool            & aSymmetry )
{
    // calculate number of samples
    uint tCount = 0;
    for(uint i=0; i < aLagrangeToBSplineMap.length(); i++)
    {
        // increment is
        // 1 sample for node id
        // 1 sample for n samples for node
        // n samples
        tCount += 2 + 2*aTMatrixFieldDesign( i, 0 );
    }

    // assign memory of output matrix
    Mat < real > tDofConnectivityFemdoc( tCount, 1 );

    // reset counter
    tCount = 0;

    if ( ! aSymmetry )
    {
        // standard case

        // loop over all Lagrange basis
        for( uint i=0; i < aLagrangeToBSplineMap.length(); ++i )
        {
            // write basis number
            tDofConnectivityFemdoc( tCount++ ) = aLagrangeToBSplineMap(i);

            // write number of samples
            tDofConnectivityFemdoc( tCount++ ) = 2*aTMatrixFieldDesign(i,0);

            // write B-Spline IDs
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = aDesignBSplineListMap.find(aIdFieldFieldDesign(i,j+1));
            }

            // write entries of T-Matrix
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = aTMatrixFieldDesign(i,j+1);
            }
        }

        // make sure that sample number is correct
        MORIS_ASSERT ( tDofConnectivityFemdoc.length() == tCount,
                "Size of tDofConnectivity is wrong" );
    }
    else
    {
        // symmetric case

        // loop over all Lagrange basis
        for( uint i=0; i < aLagrangeToBSplineMap.length(); ++i )
        {

            // calculate number of samples
            uint tNumberOfSamples = aTMatrixFieldDesign(i,0);

            // allocate memory
            Mat< uint > tIdField( tNumberOfSamples , 1 );
            Mat< real > tTMatrix ( tNumberOfSamples, 1 );

            // copy values from input
            for( uint j=0; j<tNumberOfSamples; ++j )
            {
                tIdField( j ) =  aDesignBSplineListMap.find( aIdFieldFieldDesign(i,j+1) );
                tTMatrix( j ) =  aTMatrixFieldDesign(i,j+1);
            }

            // make unique
            this-> make_design_variables_unique( tIdField, tTMatrix );

            // calculate new length
            tNumberOfSamples = tIdField.length();

            // write basis number
            tDofConnectivityFemdoc( tCount++ ) = aLagrangeToBSplineMap(i);

            // write number of samples
            tDofConnectivityFemdoc( tCount++ ) = 2*tNumberOfSamples;

            // write B-Spline IDs
            for(uint j = 0; j< tNumberOfSamples; j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = tIdField( j );
            }

            // write entries of T-Matrix
            for(uint j = 0; j< tNumberOfSamples; j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = tTMatrix( j );
            }
        }

        // reduce size of output vector
        tDofConnectivityFemdoc.resize( tCount, 1 );
    }
    // get name for output file
    std::string tOutputFile = "DesignVariables.dat";

    // pick name
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tOutputFile = "DesignVariables_new.data";
        }
        else
        {
            tOutputFile = "DesignVariables_new.data_" + std::to_string(par_rank());
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tOutputFile = "DesignVariables.data.org";
        }
        else
        {
            tOutputFile = "DesignVariables.data.org_" + std::to_string(par_rank());
        }
    }

    // write output vector to file
    save_vector_to_binary_file( tOutputFile, tDofConnectivityFemdoc );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_design_variables_for_moris_binary(
        const Mat<uint>       & aLagrangeToBSplineMap,
        const Mat<uint>       & aIdFieldFieldDesign,
        const Mat<real>       & aTMatrixFieldDesign,
        const map<uint, uint> & aDesignBSplineListMap,
        const bool            & aRefinement,
        const bool            & aSymmetry )
{
    // calculate number of samples
    uint tCount = 0;
    for(uint i=0; i < aLagrangeToBSplineMap.length(); i++)
    {
        // increment is
        // 1 sample for node id
        // 1 sample for n samples for node
        // n samples
        tCount += 2 + 2*aTMatrixFieldDesign( i, 0 );
    }

    // assign memory of output matrix
    Mat < real > tDofConnectivityFemdoc( tCount, 1 );

    // reset counter
    tCount = 0;

    if ( ! aSymmetry )
    {
        // standard case

        // loop over all Lagrange basis
        for( uint i=0; i < aLagrangeToBSplineMap.length(); ++i )
        {
            // write basis number
            tDofConnectivityFemdoc( tCount++ ) = aLagrangeToBSplineMap(i);

            // write number of samples
            tDofConnectivityFemdoc( tCount++ ) = 2*aTMatrixFieldDesign(i,0);

            // write B-Spline IDs
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = aDesignBSplineListMap.find(aIdFieldFieldDesign(i,j+1));
            }

            // write entries of T-Matrix
            for(uint j = 0; j< aIdFieldFieldDesign(i,0); j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = aTMatrixFieldDesign(i,j+1);
            }
        }

        // make sure that sample number is correct
        MORIS_ASSERT ( tDofConnectivityFemdoc.length() == tCount,
                "Size of tDofConnectivity is wrong" );
    }
    else
    {
        // symmetric case

        // loop over all Lagrange basis
        for( uint i=0; i < aLagrangeToBSplineMap.length(); ++i )
        {

            // calculate number of samples
            uint tNumberOfSamples = aTMatrixFieldDesign(i,0);

            // allocate memory
            Mat< uint > tIdField( tNumberOfSamples , 1 );
            Mat< real > tTMatrix ( tNumberOfSamples, 1 );

            // copy values from input
            for( uint j=0; j<tNumberOfSamples; ++j )
            {
                tIdField( j ) =  aDesignBSplineListMap.find( aIdFieldFieldDesign(i,j+1) );
                tTMatrix( j ) =  aTMatrixFieldDesign(i,j+1);
            }

            // make unique
            this-> make_design_variables_unique( tIdField, tTMatrix );

            // calculate new length
            tNumberOfSamples = tIdField.length();

            // write basis number
            tDofConnectivityFemdoc( tCount++ ) = aLagrangeToBSplineMap(i);

            // write number of samples
            tDofConnectivityFemdoc( tCount++ ) = 2*tNumberOfSamples;

            // write B-Spline IDs
            for(uint j = 0; j< tNumberOfSamples; j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = tIdField( j );
            }

            // write entries of T-Matrix
            for(uint j = 0; j< tNumberOfSamples; j++)
            {
                tDofConnectivityFemdoc( tCount++ ) = tTMatrix( j );
            }
        }

        // reduce size of output vector
        tDofConnectivityFemdoc.resize( tCount, 1 );
    }

    // get name for output file
    std::string tOutputFile = "DesignVariables_for_moris.data";

    // pick name
    if(  aRefinement > 0)
    {
        if( par_size() == 1)
        {
            tOutputFile = "DesignVariables_for_moris_new.data";
        }
        else
        {
            tOutputFile = "DesignVariables_for_moris_new.data_" + std::to_string(par_rank());
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tOutputFile = "DesignVariables_for_moris.data.org";
        }
        else
        {
            tOutputFile = "DesignVariables_for_moris.data.org_" + std::to_string(par_rank());
        }
    }

    // the only difference is two extra entries in the beginning

    // assign memory of output matrix
    Mat < real > tDofConnectivityMoris( tCount+2, 1 );

    tDofConnectivityMoris( 0 ) = aIdFieldFieldDesign.n_rows();
    tDofConnectivityMoris( 1 ) = aIdFieldFieldDesign.n_cols();

    // copy data from existing vector
    tDofConnectivityMoris.rows( 2, tCount+1 )
                        = tDofConnectivityFemdoc.rows( 0, tCount-1 );

    // write output vector to file
    save_vector_to_binary_file( tOutputFile, tDofConnectivityMoris );

}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::output_active_design_elements(
        Mat<uint> aElementListActiveDesign,
        bool & aRefinement)
{
    std::string tFileName;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveDesignElements_new.data";
        }
        else
        {
            tFileName = ("ActiveDesignElements_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveDesignElements.data.org";
        }
        else
        {
            tFileName = ("ActiveDesignElements_new.data.org_" + std::to_string(par_rank()));
        }
    }

    // write list of active elements into output
    save_vector_to_binary_file( tFileName, aElementListActiveDesign );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::save_active_fem_elements(
        bool & aRefinement,
        Mat<uint> & aElementListOnProc)
{
    std::string tFileName;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveFEMElements_new.data";
        }
        else
        {
            tFileName = ("ActiveFEMElements_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tFileName = ("ActiveFEMElements.data.org");
        }
        else
        {
            tFileName = ("ActiveFEMElements.data.org_" + std::to_string(par_rank()));
        }
    }

    // write list of active elements into output
    save_vector_to_binary_file( tFileName, aElementListOnProc );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::save_active_design_basis(
        bool      & aRefinement,
        Mat<uint> & aBasisActiveDesignList )
{
    std::string tFileName;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveDesignBasis_new.data";
        }
        else
        {
            tFileName = ("ActiveDesignBasis_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveDesignBasis.data.org";
        }
        else
        {
            tFileName = ("ActiveDesignBasis.data.org_" + std::to_string(par_rank()));
        }
    }

    save_vector_to_binary_file( tFileName,  aBasisActiveDesignList );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::save_active_basis(
        bool      & aRefinement,
        Mat<uint> & aBasisActiveList)
{
    std::string tFileName;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveBasis_new.data";
        }
        else
        {
            tFileName = ("ActiveBasis_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tFileName = "ActiveBasis.data.org";
        }
        else
        {
            tFileName = ("ActiveBasis.data.org_" + std::to_string(par_rank()));
        }
    }

    // write vector into file
    save_vector_to_binary_file( tFileName, aBasisActiveList );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Output::save_coordinate_list(
        const bool      & aRefinement,
        const Mat<uint> & aNodalLocaltoGlobalExist)
{
    std::string tFileName;
    if(  aRefinement == 1)
    {
        if( par_size() == 1)
        {
            tFileName = "CoordList_new.data";
        }
        else
        {
            tFileName = ("CoordList_new.data_" + std::to_string(par_rank()));
        }
    }
    else
    {
        if( par_size() == 1)
        {
            tFileName = "CoordList.data.org";
        }
        else
        {
            tFileName = ("CoordList.data.org_" + std::to_string(par_rank()));
        }
    }

    // write coordinate list into output
    save_vector_to_binary_file( tFileName, aNodalLocaltoGlobalExist );
}

// -----------------------------------------------------------------------------

void
Hierarchical_Mesh_Output::save_field_to_file(
        const std::string      & aFilePath,
        const BoostBitset      & aBasisActive,
        const map< uint, real> & aFieldMap
        )
{

    // get total number of basis functions
    uint tNumberOfBasis = aBasisActive.size();

    // get number of active basis functions
    uint tNumberOfActiveBasis = aBasisActive.count();

    // assign memory for output vector
    Mat< real > aOutput( tNumberOfActiveBasis, 1 );

    // counter for basus
    uint tCount = 0;

    // loop over all basis
    for ( uint k=0; k<tNumberOfBasis; ++k )
    {
        // check if basis is active
        if ( aBasisActive.test ( k ) )
        {
            // write value into output
            aOutput( tCount++ ) = aFieldMap.find( k );
        }
    }

    // save output into file
    save_vector_to_binary_file( aFilePath, aOutput );
}

// -----------------------------------------------------------------------------

void
Hierarchical_Mesh_Output::make_design_variables_unique(
        Mat< uint > & aIdField,
        Mat< real > & aTMatrix )
{
    // make ID field unique
    Mat< uint > tIdField = unique( aIdField );

    // create small map
    map< uint, uint > tMap;

    // counter for Map
    uint tCount = 0;
    for ( uint k=0; k< tIdField.length(); ++k )
    {
        tMap[ tIdField( k ) ] = tCount++;
    }

    // create output T-Matrix
    Mat< real > tTMatrix( tCount, 1, 0);

    // assemble new T-Matrix
    for ( uint k=0; k< aIdField.length(); ++k )
    {
        tTMatrix( tMap.find( aIdField( k ) ) ) +=  aTMatrix( k ) ;
    }

    // overwrite output variables
    aIdField = tIdField;
    aTMatrix = tTMatrix;
}

// -----------------------------------------------------------------------------
