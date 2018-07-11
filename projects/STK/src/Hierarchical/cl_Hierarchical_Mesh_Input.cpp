/*
 * cl_Hierarchical_Mesh_Input.cpp
 *
 *  Created on: Jan 8, 2018
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Input.hpp" // STK/src/Hierarchical
using namespace moris;

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_active_fem_elements(
        Mat<uint> & aElemLocaltoGlobalLastStepFEM,
        uint & aLevel,
        const uint & aDim,
        const Mat<uint> & aNumElements,
        BoostBitset & aElementActiveLastStep)
{

    // determine filename
    std::string tFileName;
    if( par_size() == 1)
    {
        tFileName = ("ActiveFEMElements.data");
    }
    else
    {
        tFileName = ("ActiveFEMElements.data_" + std::to_string(par_rank()));
    }

    // read binary file
    load_vector_from_binary_file( tFileName, aElemLocaltoGlobalLastStepFEM );

    // calculate maximum level of mesh
    aLevel = mBaseElement.give_element_level(
            aElemLocaltoGlobalLastStepFEM.max(),
            aDim,
            aNumElements);

    // calculate size of bitset
    uint tNumberBasis = mBaseElement.give_number_of_elements(
            aLevel,
            aDim,
            aNumElements);
    // create bitset from vector
    aElementActiveLastStep.resize( tNumberBasis );
    aElementActiveLastStep.reset();
    uint tNumberOfElements = aElemLocaltoGlobalLastStepFEM.length();
    for( uint i = 0; i < tNumberOfElements; ++i )
    {
        aElementActiveLastStep.set( aElemLocaltoGlobalLastStepFEM( i ) );
    }

    std::fprintf(stdout, "Number of active FEM elements read: %u\n", tNumberOfElements );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_active_design_elements(
        const uint      & aDim,
        const Mat<uint> & aNumElements,
        uint            & aLevel,
        Mat<uint>       & aElemLocaltoGlobalLastStepDesign,
        BoostBitset     & aElementActiveDesignLastStep)
{

    // determine filename
    std::string tFileName;
    if( par_size() == 1)
    {
        tFileName = ("ActiveDesignElements.data");
    }
    else
    {
        tFileName = ("ActiveDesignElements.data_" + std::to_string(par_rank()));
    }

    // load binary file
    load_vector_from_binary_file( tFileName, aElemLocaltoGlobalLastStepDesign );

    // calculate maximum refinement level
    aLevel = mBaseElement.give_element_level(
            aElemLocaltoGlobalLastStepDesign.max(),
            aDim,
            aNumElements);

    // calculate size of bitset
    uint tNumberElements = mBaseElement.give_number_of_elements(
            aLevel,
            aDim,
            aNumElements);

    // create bitset form vector
    aElementActiveDesignLastStep.resize( tNumberElements );
    aElementActiveDesignLastStep.reset();

    uint tNumberOfElements = aElemLocaltoGlobalLastStepDesign.length();
    for( uint i = 0; i < tNumberOfElements; ++i )
    {
        aElementActiveDesignLastStep.set( aElemLocaltoGlobalLastStepDesign( i ) );
    }

    std::fprintf(stdout, "Number of active design elements read: %u\n", tNumberOfElements );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_active_design_basis(
        const uint      & aDim,
        const uint      & aPolynomialDesign,
        const Mat<uint> & aNumElements,
        Mat<uint>       & aBasisActiveDesignListLastStep,
        BoostBitset     & aBasisActiveDesign )
{
    // determine filename
    std::string tFileName;
    if( par_size() == 1)
    {
        tFileName = ("ActiveDesignBasis.data");
    }
    else
    {
        tFileName = ("ActiveDesignBasis.data_" + std::to_string(par_rank()));
    }

    // load binary file
    load_vector_from_binary_file( tFileName, aBasisActiveDesignListLastStep );

    // calculate maximum design level
    uint tLevel = mBasis.give_basis_level(
            aBasisActiveDesignListLastStep.max(),
            aDim,
            aPolynomialDesign,
            aNumElements) + 1;

    // calculate size of bitset
    uint tNumberBasis = mBasis.give_number_of_basis(
            tLevel,
            aDim,
            aPolynomialDesign,aNumElements );

    // create bitset from vector
    aBasisActiveDesign.resize( tNumberBasis );
    aBasisActiveDesign.reset();

    uint tNumberOfBasis = aBasisActiveDesignListLastStep.length();
    for(uint i = 0; i < tNumberOfBasis; ++i)
    {
        aBasisActiveDesign.set( aBasisActiveDesignListLastStep( i ) );
    }

    std::fprintf(stdout, "Number of active design basis read: %u\n", tNumberOfBasis );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_active_basis(
        const uint      & aDim,
        const uint      & aPolynomialDesign,
        const Mat<uint> & aNumElements,
        Mat<uint>       & aBasisActiveListLastStep,
        BoostBitset     & aBasisActive)
{

    // determine filename
    std::string tFileName;
    if( par_size() == 1)
    {
        tFileName = ("ActiveBasis.data");
    }
    else
    {
        tFileName = ("ActiveBasis.data_" + std::to_string(par_rank()));
    }

    // load binary file
    load_vector_from_binary_file( tFileName, aBasisActiveListLastStep );
    // calculate max level of basis
    uint tLevel = mBasis.give_basis_level(
            aBasisActiveListLastStep.max(),
            aDim,
            aPolynomialDesign,
            aNumElements) + 1;

    // calculate sizs of bitset
    uint tNumberBasis = mBasis.give_number_of_basis(
            tLevel,
            aDim,
            aPolynomialDesign,
            aNumElements);

    // create bitset from vector
    aBasisActive.resize( tNumberBasis );
    aBasisActive.reset();
    uint tNumberOfBasis = aBasisActiveListLastStep.length();
    for( uint i = 0; i < tNumberOfBasis; ++i )
    {
        aBasisActive.set( aBasisActiveListLastStep( i ) );
    }

    std::fprintf(stdout, "Number of active FEM basis read: %u\n", tNumberOfBasis );
}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_design_variables_ascii(
        Mat<uint> & aIdFieldDesignNodalField,
        Mat<real> & aTMatrixDesignNodalField)
{
    //Read active elements from file
    uint tNcolsIdField;
    uint tNrowsIdField;
    uint tDummyUINTValLoop;
    uint tDummyUINTVal;
    real tDummyREALVal;
    std::ifstream infileee("DesignVariables_for_moris.data");
    MORIS_ASSERT( infileee, "Cannot open DesignVariables_for_moris.data\n" );
    infileee >> tNrowsIdField;
    infileee >> tNcolsIdField;

    aIdFieldDesignNodalField.set_size(tNrowsIdField,tNcolsIdField,0);
    aTMatrixDesignNodalField.set_size(tNrowsIdField,tNcolsIdField,0);
    for(uint i = 0; i < tNrowsIdField; i++)
    {
        infileee >> tDummyUINTValLoop;
        infileee >> tDummyUINTValLoop;
        tDummyUINTValLoop /=2;
        aIdFieldDesignNodalField(i,0) = tDummyUINTValLoop;
        for(uint j = 0; j < tDummyUINTValLoop; j++)
        {
            infileee >> tDummyUINTVal;
            aIdFieldDesignNodalField(i,j+1) = tDummyUINTVal;
        }
        aTMatrixDesignNodalField(i,0) = tDummyUINTValLoop;
        for(uint j = 0; j < tDummyUINTValLoop; j++)
        {
            infileee >> tDummyREALVal;
            aTMatrixDesignNodalField(i,j+1) = tDummyREALVal;
        }
    }

    std::fprintf(stdout, "Number of design variables read:  %u\n", tNrowsIdField );

}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_design_variables_binary(
        Mat<uint> & aIdFieldDesignNodalField,
        Mat<real> & aTMatrixDesignNodalField)
{
    //Read active elements from file
    uint tNumberOfSamples;

    std::string tFilePath = "DesignVariables_for_moris.data";

    // the matrix containing the data in the file
    Mat< real > aData;

    load_vector_from_binary_file( tFilePath, aData );

    // reset counter
    uint tCount = 0;

    // read rows and cols form file
    uint tNrowsIdField = ( uint ) aData( tCount++ );
    uint tNcolsIdField = ( uint ) aData( tCount++ );

    aIdFieldDesignNodalField.set_size( tNrowsIdField, tNcolsIdField, 0 );
    aTMatrixDesignNodalField.set_size( tNrowsIdField, tNcolsIdField, 0 );

    for( uint i = 0; i < tNrowsIdField; ++i )
    {
        // skip Lagrange ID for FEMDOC ( using internal MORIS ID instead )
        ++tCount;

        // read number of samples
        tNumberOfSamples = aData( tCount++ );

        // divide by two
        tNumberOfSamples /= 2;

        // write number of samples in first column of ID field
        aIdFieldDesignNodalField(i,0) = tNumberOfSamples;

        // read ID field
        for(uint j = 0; j < tNumberOfSamples; ++j )
        {
            aIdFieldDesignNodalField(i,j+1) = ( uint ) aData( tCount++ );
        }

        // write number of samples in first column if T-Matrix
        aTMatrixDesignNodalField(i,0) = tNumberOfSamples;
        for(uint j = 0; j < tNumberOfSamples; ++j )
        {
            aTMatrixDesignNodalField(i,j+1) =  aData( tCount++ );;
        }
    }

    // make sure that all data is used
    MORIS_ASSERT( tCount == aData.length(), "Something went wrong while reading from file DesignVariables_for_moris.data." );

    std::fprintf(stdout, "Number of design variables read: %u\n", tNrowsIdField );

}

// -----------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_coordinate_list(
        Mat<uint> & aCoordinateList)
{
    std::string tFileName;
    if( par_size() == 1)
    {
        tFileName = ("CoordList.data");
    }
    else
    {
        tFileName = ("CoordList.data_" + std::to_string(par_rank()));
    }

    load_vector_from_binary_file( tFileName,  aCoordinateList );

    std::fprintf(stdout, "number of entries read from CoordList.data: %u\n",
           ( uint ) aCoordinateList.length() );
}

// -----------------------------------------------------------------------------
// FIXME: This function is obsolete. Remove it.
void Hierarchical_Mesh_Input::read_sdf_file(
        std::string const & aSDFFileName,
        Mat<real> & aSDFNodalField,
        Mat<uint> & aCoordinateList)
{
    MORIS_ASSERT( isempty(aCoordinateList)==0,"Need to provide a Coordinate list of the mesh to put the SDF file solutions on the right positions");
    char *cstr = new char[aSDFFileName.length() + 1];
    strcpy(cstr, aSDFFileName.c_str());
    // do stuff
    FILE* fid=std::fopen(cstr,"r");
    delete [] cstr;
    MORIS_ASSERT( fid != NULL, "Cannot open SDF File\n" );
    fseek(fid, 0L, SEEK_END);
    uint tNumCoordinates;
    tNumCoordinates = ftell(fid)/sizeof(uint);
    fseek(fid, 0L, SEEK_SET);
    uint * buffer = (uint*) alloca(tNumCoordinates*sizeof(uint));
    std::fread(buffer,sizeof(uint),tNumCoordinates,fid);
    std::fclose(fid);
    aSDFNodalField.set_size(tNumCoordinates,1,0);
    MORIS_ASSERT( aCoordinateList.length() == tNumCoordinates,"Number of Coordinates are different from those of the SDF file");
    for(uint i = 0; i < tNumCoordinates; i++)
        aCoordinateList(aCoordinateList(i)) = buffer[i];
}

//-------------------------------------------------------------------------------
// Read design variables defined on design basis of last step and
// create level set, density, and nodal ADV fields on last active basis mesh
//-------------------------------------------------------------------------------

void Hierarchical_Mesh_Input::read_absdesvariables_file(
        const uint                      & aDim,
        const uint                      & aPolynomialDesign,
        const Mat<uint>                 & aNumElements,
        const Mat<uint>                 & aNumBasis,
        const Mat<uint>                 & aBasisActiveDesignListLastStep,
        map<uint, real>                 & aNodalADV,
        map<uint, real>                 & aNodalPDV,
        map<uint, real>                 & aNodalLvlsetField,
        const Cell< map< uint, real > > & aNodalSDFs,
        BoostBitset                     & aNodalFieldExists,
        Mat<uint>                       & tIdFieldDesignNodalField,
        Mat<real>                       & tTMatrixDesignNodalField,
        real                            & tSimpExp,
        real                            & sElementEdgeWidth,
        real                            & sInitialElementEdgeWidth,
        real                            & sLSscale,
        real                            & sLSthresh,
        real                            & sDensityShift,
        real                            & sPrjbeta ,
        const Mat< real >               & aSDFAlpha,
        const Mat< real >               & aSDFBeta,
        const bool                      & aExpectBinaryInput,
        const bool                      & aUseSymmetry,
        const uint                      & aSymmetryPlane,
        const uint                      & aSymmetryIndex)
{
    // field containing design variables
    Mat<real> tADVField;

    // load ADV Field from binary path
    load_vector_from_binary_file( "AbsDesVariables.dat" , tADVField );

    if ( aUseSymmetry )
    {
        this->expand_symmetric_field(
                tADVField,
                aBasisActiveDesignListLastStep,
                aDim,
                aPolynomialDesign,
                aSymmetryPlane,
                aSymmetryIndex,
                aNumElements,
                aNumBasis );
    }

    uint numvar = tADVField.length();

    // check that number of variables is equal to the number of active design basis of last mesh
    MORIS_ASSERT( numvar == aBasisActiveDesignListLastStep.length(), "AbsDesVariables.dat and ActiveDesignBasis.data have not the same number of entries ." );

    // Projection function variables
    real sPrjeta  = 0.0001;
    real tBulkMin = 1e-6;
    real tProBulk = 1.0 ;
    real tLSscale = sLSscale * sInitialElementEdgeWidth;

    // Allocate vectors: fields at active basis nodes
    Mat<real> tNodalLevelSet(numvar,1);
    Mat<real> tNodalDensity(numvar,1);
    // Mat<real> tADVField(numvar,1);

    // Read the nodal ID-field and T-Matrix from a file
    // tIdFieldDesignNodalField: node ID of nodes in last Lagrange mesh
    // tTMatrixDesignNodalField: interpolation weights for node in last Lagrange mesh
    if ( aExpectBinaryInput )
    {
        this->read_design_variables_binary( tIdFieldDesignNodalField, tTMatrixDesignNodalField );
    }
    else
    {
        this->read_design_variables_ascii( tIdFieldDesignNodalField, tTMatrixDesignNodalField );
    }

    if ( ! aUseSymmetry )
    {
        MORIS_ASSERT( (tIdFieldDesignNodalField.max()+1) == numvar, "AbsDesVariables.dat has not the same size as the Designvariables in Designvariables_for_moris.data");
    }

    // Fixme: the calculations below are not consistent with what is done in femdoc; this needs to be a loop over all nodes in the mesh; not a loop over all variables
    //Calculate from the AbsDesVariables file a user specific nodal field
    for (uint i = 0; i < numvar; i++)
    {
        //real avePara = 0.0;
        // FIXME: find out why these lines do not work
        //        for(uint j = 0; j < tIdFieldDesignNodalField(i,0); j++)
        //        {
        //            avePara += buffer[tIdFieldDesignNodalField(i,j+1)] * tTMatrixDesignNodalField(i,j+1);    //(buffer[i]-sLSshift)/(1-sLSshift);
        //        }
        //avePara      = buffer[i];
        //tADVField(i) = avePara;

        real lsValue;
        real density;

        real prjv=std::tanh(sPrjbeta*sPrjeta);
        real prjq=std::tanh(sPrjbeta*sPrjeta)+std::tanh(sPrjbeta*(1-sPrjeta));
        real valp=(std::tanh(sPrjbeta*(tADVField(i)-sPrjeta))+prjv)/prjq;

        uint tNumberOfSDFs =  aNodalSDFs.size();
        if ( tNumberOfSDFs > 0 )
        {
            Mat < real > tSDF( tNumberOfSDFs, 1);

            // loop over all SDFs
            for( uint k=0; k<tNumberOfSDFs; ++k )
            {
                tSDF( k ) = aSDFAlpha( k ) * aNodalSDFs( k ).find( aBasisActiveDesignListLastStep( i ) )  + aSDFBeta( k );
            }

            // point within box or bolt and in vicinity of their boundaries
            if ( ( tSDF(1) < 0 && tSDF(1) > -sElementEdgeWidth ) || ( tSDF(2) < 0 && tSDF(2) > -sElementEdgeWidth ) )
            {
                lsValue = std::min( tSDF(1), tSDF(2) );
                density = 1.0;
            }

            // point within box or bolt but not in vicinity of their boundaries
            else if ( tSDF(1) <= -sElementEdgeWidth || tSDF(2) <= -sElementEdgeWidth )
            {
                lsValue = std::max( -tSDF(1), -tSDF(2) );
                density = 0.0;
            }

            // point outside design domain
            else if ( tSDF(0) > 0 )
            {
                lsValue = tSDF(0);
                density = 0.0;
            }

            // point within design domain
            else
            {
                lsValue =  tLSscale*(sLSthresh-valp);
                density = -lsValue/tLSscale/(1-sLSthresh);

                if ( density < 0 ) density = 0;
                if ( density > 1 ) density = 1;

                density = sDensityShift + (1-sDensityShift) * density;
            }
        }
        else
        {
            lsValue = tLSscale*(sLSthresh-valp);
            density = -lsValue/tLSscale/(1-sLSthresh);

            if ( density < 0 ) density = 0;
            if ( density > 1 ) density = 1;

            density = sDensityShift + (1-sDensityShift) * density;
        }

        // save level set and density values
        tNodalLevelSet(i) = lsValue;
        tNodalDensity(i)  = tBulkMin + (tProBulk-tBulkMin)*std::pow(density,tSimpExp);
    }

    // reset vectors defining fields with respect to ID of active design basis of last mesh
    aNodalADV.clear();
    aNodalPDV.clear();
    aNodalLvlsetField.clear();

    aNodalFieldExists.resize(aBasisActiveDesignListLastStep.max() + 1);
    aNodalFieldExists.reset();
    for(uint i = 0; i < aBasisActiveDesignListLastStep.length(); i++)
    {
        // mark design basis active in bitset of all design basis of last step
        aNodalFieldExists.set(aBasisActiveDesignListLastStep(i));

        // populate maps; key: design basis node ID
        aNodalADV[aBasisActiveDesignListLastStep(i)]    = tADVField(i);
        aNodalPDV[aBasisActiveDesignListLastStep(i)]    = tNodalDensity(i);
        aNodalLvlsetField[aBasisActiveDesignListLastStep(i)] = tNodalLevelSet(i);
    }

    std::fprintf(stdout, "number of variables read from AbsDesVariables.dat: %u\n",
            ( uint ) aBasisActiveDesignListLastStep.length() );
}

// -----------------------------------------------------------------------------

void
Hierarchical_Mesh_Input::read_field_from_file(
        const std::string & aFilePath,
        const BoostBitset & aBasisActive,
        const bool        & aUseSymmetry,
        const Mat< uint > & aBasisActiveDesignListLastStep,
        const uint        & aModelDim,
        const uint        & aPolynomial,
        const uint        & aSymmetryPlane,
        const uint        & aSymmetryIndex,
        const Mat<uint>   & aNumberOfElementsPerDirection,
        const Mat<uint>   & aNumberOfBasisPerDirection,
        map< uint, real>  & aFieldMap )
{

    // get total number of basis functions
    uint tNumberOfBasis = aBasisActive.size();

    // assign memory for input
    Mat< real > tInput;

    // read vector from file
    load_vector_from_binary_file( aFilePath, tInput );

    if ( aUseSymmetry )
    {
        this->expand_symmetric_field(
                        tInput,
                        aBasisActiveDesignListLastStep,
                        aModelDim,
                        aPolynomial,
                        aSymmetryPlane,
                        aSymmetryIndex,
                        aNumberOfElementsPerDirection,
                        aNumberOfBasisPerDirection );
    }

    // make sure that input vector has the right length
    MORIS_ASSERT( tInput.length() == aBasisActive.count(),
               "Length of field file to read from file does not match number of active basis.");

    // counter for basis
    uint tCount = 0;

    // reset output map
    aFieldMap.clear();

    // loop over all basis
    for ( uint k=0; k<tNumberOfBasis; ++k )
    {
        // check if basis is active
        if ( aBasisActive.test ( k ) )
        {
            // write value into field
            aFieldMap[ k ] = tInput( tCount++ );
        }
    }
}

// -----------------------------------------------------------------------------

void
Hierarchical_Mesh_Input::expand_symmetric_field(
        Mat< real >       & aField,
        const Mat< uint > & aBasisActiveDesignListLastStep,
        const uint        & aModelDim,
        const uint        & aPolynomial,
        const uint        & aSymmetryPlane,
        const uint        & aSymmetryIndex,
        const Mat<uint>   & aNumberOfElementsPerDirection,
        const Mat<uint>   & aNumberOfBasisPerDirection )
{

    // make sure that we are working in serial
    MORIS_ASSERT( par_size() == 1, "Hierarchical_Mesh_Input::expand_symmetric_field not tested for parallel" );

    // count active basis
    uint tNumberOfActiveBasis = aBasisActiveDesignListLastStep.length();
    // counter
    uint tCount = 0;

    // step 1: generate order in which data in field is stored
    Mat< uint > tSymmetricMaster( tNumberOfActiveBasis, 1 );
    Mat< uint > tSymmetricSlave( tNumberOfActiveBasis, 1 );
    // loop over all active basis
    for ( uint k=0; k< tNumberOfActiveBasis; ++k )
    {
        // get ID of basis
        uint tThisBasis = aBasisActiveDesignListLastStep( k );

        // check if basis is symmetry master
        if ( mBasis.basis_is_symmetry_master(
                tThisBasis,
                aModelDim,
                aPolynomial,
                aSymmetryPlane,
                aSymmetryIndex,
                aNumberOfElementsPerDirection ) )
        {
            // add basis to master list
            tSymmetricMaster( tCount ) = tThisBasis;

            // calculate slave ( note, in parallel, we must test here if this node is shared by current proc )
            tSymmetricSlave( tCount ) = mBasis.give_symmetry_slave_of_basis(
                    tThisBasis,
                    aModelDim,
                    aPolynomial,
                    aSymmetryPlane,
                    aSymmetryIndex,
                    aNumberOfElementsPerDirection,
                    aNumberOfBasisPerDirection ) ;

            // increment counter
            ++tCount;
        }
    }

    // resize arrays
    tSymmetricMaster.resize( tCount, 1);
    tSymmetricSlave.resize( tCount, 1);

    // save for debugging
    //save_vector_to_binary_file( "symmetric_master.mat", tSymmetricMaster );

    // Step 2: create map for active basis from last step

    map< uint, uint > tMap;

    // loop over all basis
    for ( uint k=0; k< tNumberOfActiveBasis; ++k )
    {
        tMap[ aBasisActiveDesignListLastStep( k ) ] = k;
    }
    // Step 3: expand vector to all basis
    Mat< real > tField( tNumberOfActiveBasis, 1);

    // Boost bitset for debugging
    BoostBitset tChecklist( tNumberOfActiveBasis );

    // loop over all symmetric basis basis
    for ( uint k=0; k<tCount; ++k )
    {

        // get position of master
        uint tMaster = tMap.find( tSymmetricMaster( k ) );

        // get position of slave
        uint tSlave = tMap.find( tSymmetricSlave( k ) );

        // copy input data into new field
        tField( tMaster ) = aField( k );
        tField( tSlave )  = aField( k );

        // check master and slave
        tChecklist.set( tMaster );
        tChecklist.set( tSlave  );
    }

    // make sure that we got all of them
    MORIS_ASSERT( tChecklist.count() == tNumberOfActiveBasis,
                  "something went wrong win Hierarchical_Mesh_Input::expand_symmetric_field" );

    // return output  vector
    aField = tField;

}
