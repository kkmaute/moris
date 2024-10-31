/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Meshgen.cpp
 *
 */

#include "cl_Library_IO_Meshgen.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_stringify_matrix.hpp"
#include "cl_XML_Parser.hpp"
#include "enums.hpp"
#include "parameters.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Meshgen::Library_IO_Meshgen()
            : Library_IO()    // initialize base class data as usual
    {
        // list of supported parameter list types
        mSupportedParamListTypes = { Module_Type::HMR, Module_Type::XTK, Module_Type::GEN };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Meshgen::~Library_IO_Meshgen()
    {
        // do nothing extra
    }

    // -----------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_parameters_from_xml()
    {
        // check that an XML file has been initialized
        MORIS_ERROR( mXmlParserIsInitialized, "Library_IO::load_parameters_from_xml() - No XML file has been loaded." );

        // declare root of the input file and the parameter lists
        std::string tXmlRoot = "MeshGenerationParameterList";
        std::string tGenRoot = "Geometries";
        std::string tHmrRoot = "BackgroundMeshes";
        std::string tXtkRoot = "ForegroundMesh";

        // get the paths for the various parameter lists
        std::string tGenPath = tXmlRoot + "." + tGenRoot;
        std::string tHmrPath = tXmlRoot + "." + tHmrRoot;
        std::string tXtkPath = tXmlRoot + "." + tXtkRoot;

        // check if all necessary modules have been declared
        size_t tGenCount = mXmlReader->count_keys_in_subtree( tXmlRoot, "Geometries" );
        size_t tHmrCount = mXmlReader->count_keys_in_subtree( tXmlRoot, "BackgroundMeshes" );
        size_t tXtkCount = mXmlReader->count_keys_in_subtree( tXmlRoot, "ForegroundMesh" );

        MORIS_ERROR( tGenCount == 1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Parameters for 'Geometries' either missing or declared multiple times." );
        MORIS_ERROR( tHmrCount == 1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Parameters for 'BackgroundMeshes' either missing or declared multiple times." );
        MORIS_ERROR( tXtkCount == 1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Parameters for 'ForegroundMesh' either missing or declared multiple times." );

        // load the parameters
        this->load_HMR_parameters_from_xml( tHmrPath, tXtkPath );
        this->load_XTK_parameters_from_xml( tXtkPath, tHmrPath );
        this->load_GEN_parameters_from_xml( tGenPath, tHmrPath, tXtkPath );

    }    // end function: Library_IO_Meshgen::load_parameters_from_xml()

    // -----------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_HMR_parameters_from_xml(
            std::string const & aHmrPath,
            std::string const & aXtkPath )
    {
        // quick access to the parameter list
        Parameter_List& tHmrParamList = mParameterLists( (uint)( Module_Type::HMR ) )( 0 )( 0 );

        // ------------------------------
        // Base grid

        // get the base grid parameters
        std::string tBaseGridPath = aHmrPath + ".BaseGrid";
        std::string tBaseGridSize;
        std::string tDomainDimensions;
        std::string tBaseGridOrigin;
        mXmlReader->get( tBaseGridPath + ".Size", tBaseGridSize, std::string( "" ) );
        mXmlReader->get( tBaseGridPath + ".Dimensions", tDomainDimensions, std::string( "" ) );
        mXmlReader->get( tBaseGridPath + ".Origin", tBaseGridOrigin, std::string( "" ) );

        // check the number of spatial dimensions
        Matrix< DDUMat > tBaseGridMat;
        moris::string_to_matrix( tBaseGridSize, tBaseGridMat );
        mNumSpatialDims = tBaseGridMat.numel();

        Matrix< DDRMat > tDomainDimsMat;
        moris::string_to_matrix( tDomainDimensions, tDomainDimsMat );
        uint tNumDimsDimensions = tDomainDimsMat.numel();

        Matrix< DDRMat > tDomainOriginMat;
        moris::string_to_matrix( tBaseGridOrigin, tDomainOriginMat );
        uint tNumDimsOrigin = tDomainOriginMat.numel();

        // perform some checks on the user inputs
        MORIS_ERROR( mNumSpatialDims == 2 || mNumSpatialDims == 3,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid specified must have 2 or 3 dimensions." );
        MORIS_ERROR( mNumSpatialDims == tNumDimsDimensions,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid 'Dimensions' inconsistent with 'Size'." );
        MORIS_ERROR( mNumSpatialDims == tNumDimsOrigin,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid 'Origin' inconsistent with 'Size'." );

        // generate the correct side set names for this number of dimensions
        std::string tDomainSideSets = ( mNumSpatialDims == 3 ) ? "1,2,3,4,5,6" : "1,2,3,4";

        // set the parameters in the parameter lists
        tHmrParamList.set( "number_of_elements_per_dimension", tBaseGridSize );
        tHmrParamList.set( "domain_dimensions", tDomainDimensions );
        tHmrParamList.set( "domain_offset", tBaseGridOrigin );
        tHmrParamList.set( "domain_sidesets", tDomainSideSets );

        // ------------------------------
        // grids and refinements

        // path to the grids
        std::string tGridsPath    = aHmrPath + ".MeshGrids";
        std::string tGridNodeName = "MeshGrid";

        // see how many grids are specified in the input file
        uint tNumGridsSpecified = mXmlReader->count_keys_in_subtree( tGridsPath, tGridNodeName );

        // initialize a map that stores in what order the grids are defined
        std::map< moris_index, uint > tGridIndexMap;    // map: grid index | position in list
        Vector< moris_index >         tGridIndices( tNumGridsSpecified, 0 );
        Vector< moris_index >         tInitialRefinements( tNumGridsSpecified, 0 );
        Vector< moris_index >         tBoundaryRefinements( tNumGridsSpecified, 0 );
        moris_index                   tMaxGridIndex = -1;

        // go through grids specified and get their data
        for ( uint iMeshGrid = 0; iMeshGrid < tNumGridsSpecified; iMeshGrid++ )
        {
            // store the parameters for this mesh grid in the buffer
            mXmlReader->copy_subtree_into_buffer( tGridsPath, tGridNodeName, iMeshGrid );

            // get the index of this mesh grid
            moris_index tGridIndex = mXmlReader->get_index_of_buffer_subtree();

            // check that the index has been set
            MORIS_ERROR( tGridIndex != -1,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "The %i-th mesh grid specified has no index attribute. Please set an 'ind' attribute in the input file.",
                    iMeshGrid );

            // determine max index
            tMaxGridIndex = std::max( tMaxGridIndex, tGridIndex );
            MORIS_ERROR( tGridIndex < 6,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Maximum grid index exceeded. Don't specify grid indices above 5." );

            // check that a grid index is defined at most once
            MORIS_ERROR( tGridIndexMap.find( tGridIndex ) == tGridIndexMap.end(),
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Grid index %i defined twice. Please use unique indices.",
                    tGridIndex );

            // add to map
            tGridIndexMap[ tGridIndex ] = iMeshGrid;
            tGridIndices( iMeshGrid )   = tGridIndex;

            // get and store the number of refinements
            moris_index& tNumInitRefines = tInitialRefinements( iMeshGrid );
            mXmlReader->get_from_buffer( "InitialRefinements", tNumInitRefines, moris_index( 0 ) );
            moris_index& tNumBdRefines = tBoundaryRefinements( iMeshGrid );
            mXmlReader->get_from_buffer( "InterfaceRefinements", tNumBdRefines, moris_index( 0 ) );
        }

        // make sure there aren't any grid indices left unidentified
        MORIS_ERROR( tMaxGridIndex == (moris_index)(tNumGridsSpecified)-1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Some grid indices have not been defined. The maximum grid index used is %i, but only %i grids are defined.",
                tMaxGridIndex,
                tNumGridsSpecified );

        // ------------------------------
        // B-spline meshes

        // path to the meshes
        std::string tBspMeshesPath   = aHmrPath + ".BsplineMeshes";
        std::string tBspMeshNodeName = "BsplineMesh";

        // see how many B-spline meshes are specified in the input file
        uint tNumBspMeshesSpecified = mXmlReader->count_keys_in_subtree( tBspMeshesPath, tBspMeshNodeName );

        // initialize a map that stores in what order the grids are defined
        std::map< moris_index, uint > tBspMeshIndexMap;    // map: B-spline mesh index | position in list
        Vector< moris_index >         tBspMeshIndices( tNumBspMeshesSpecified, -1 );
        Vector< moris_index >         tGridsForMeshes( tNumBspMeshesSpecified, -1 );
        Vector< moris_index >         tPolyOrders( tNumBspMeshesSpecified, -1 );
        moris_index                   tMaxBspMeshIndex = -1;
        moris_index                   tMaxPolyOrder    = -1;

        // check maximum refinement levels for B-spline meshes
        moris_index tMaxUniformRefineLvlOnBspMeshes = 0;
        moris_index tMaxTotalRefineLvlOnBspMeshes   = 0;

        // go through grids specified and get their data
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshesSpecified; iBspMesh++ )
        {
            // store the parameters for this mesh grid in the buffer
            mXmlReader->copy_subtree_into_buffer( tBspMeshesPath, tBspMeshNodeName, iBspMesh );

            // get the index of this mesh grid
            moris_index tBspMeshIndex = mXmlReader->get_index_of_buffer_subtree();

            // check that the index has been set
            MORIS_ERROR( tBspMeshIndex != -1,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "The %i-th B-spline mesh specified has no index attribute. Please set an 'ind' attribute in the input file.",
                    iBspMesh );

            // determine max index
            tMaxBspMeshIndex = std::max( tMaxBspMeshIndex, tBspMeshIndex );
            MORIS_ERROR( tBspMeshIndex < 6,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Maximum B-spline mesh index exceeded. Don't specify B-spline mesh indices above 5." );

            // check that a grid index is defined at most once
            MORIS_ERROR( tBspMeshIndexMap.find( tBspMeshIndex ) == tBspMeshIndexMap.end(),
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "B-spline mesh index %i is defined twice. Please use unique indices.",
                    tBspMeshIndex );

            // add to map
            tBspMeshIndexMap[ tBspMeshIndex ] = iBspMesh;
            tBspMeshIndices( iBspMesh )       = tBspMeshIndex;

            // get and store the polynomial order
            moris_index& tPolyOrder = tPolyOrders( iBspMesh );
            mXmlReader->get_from_buffer( "PolynomialOrder", tPolyOrder, moris_index( 1 ) );
            MORIS_ERROR( tPolyOrder > 0 && tPolyOrder < 4,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Only polynomial orders 1, 2, or 3 supported." );

            // determine max index
            tMaxPolyOrder = std::max( tMaxPolyOrder, tPolyOrder );

            // get and store the corresponding grid index
            moris_index& tGridIndexUsed = tGridsForMeshes( iBspMesh );
            mXmlReader->get_from_buffer( "MeshGridIndex", tGridIndexUsed, moris_index( 0 ) );
            MORIS_ERROR( tGridIndexUsed > -1 && tGridIndexUsed < tMaxGridIndex + 1,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Trying to use mesh grid index %i for B-spline mesh. But only indices 0 to %i are defined.",
                    tGridIndexUsed,
                    tMaxGridIndex );

            // get the refinement levels used
            uint        tGridPosInMap       = tGridIndexMap[ tGridIndexUsed ];
            moris_index tInitialRefineLvl   = tInitialRefinements( tGridPosInMap );
            moris_index tBdRefines          = tBoundaryRefinements( tGridPosInMap );
            tMaxUniformRefineLvlOnBspMeshes = std::max( tMaxUniformRefineLvlOnBspMeshes, tInitialRefineLvl );
            tMaxTotalRefineLvlOnBspMeshes   = std::max( tMaxTotalRefineLvlOnBspMeshes, tInitialRefineLvl + tBdRefines );
        }

        // make sure there aren't any mesh indices left unidentified
        MORIS_ERROR( tMaxBspMeshIndex == (moris_index)(tNumBspMeshesSpecified)-1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Some B-spline mesh indices have not been defined. The maximum mesh index used is %i, but only %i meshes are defined.",
                tMaxBspMeshIndex,
                tNumBspMeshesSpecified );

        // ------------------------------
        // Lagrange mesh

        // get which grid is used for decomposition
        moris_index tGridIndexForDecomp = 0;
        mXmlReader->get( aXtkPath + ".DecompositionGrid", tGridIndexForDecomp, 0 );
        MORIS_ERROR( tGridIndexForDecomp > -1 && tGridIndexForDecomp < tMaxGridIndex + 1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Trying to use mesh grid index %i for decomposition. But only grid indices 0 to %i are defined.",
                tGridIndexForDecomp,
                tMaxGridIndex );

        // make sure the grid specified for decomposition actually exists in the background mesh parameter list
        auto tDecompGridIter = tGridIndexMap.find( tGridIndexForDecomp );
        MORIS_ERROR( tDecompGridIter != tGridIndexMap.end(),
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Grid index for decomposition not found in map. Something went wrong." );

        // make sure the Lagrange mesh is the most refined one (if not create a new grid for it later)
        uint        tDecompGridPosInList      = tDecompGridIter->second;
        moris_index tUniformRefinesOnFg       = tInitialRefinements( tDecompGridPosInList );
        moris_index tTotalRefinesOnFg         = tUniformRefinesOnFg + tBoundaryRefinements( tDecompGridPosInList );
        bool        tCreateMoreRefinedLagMesh = ( tUniformRefinesOnFg < tMaxUniformRefineLvlOnBspMeshes ) || ( tTotalRefinesOnFg < tMaxTotalRefineLvlOnBspMeshes );

        // add data for additional grid for decomposition if it needs to
        if ( tCreateMoreRefinedLagMesh )
        {
            // get new grid's index and set the decomp grid to it
            tMaxGridIndex++;
            tGridIndexForDecomp = tMaxGridIndex;

            // add grid index to map
            tGridIndexMap[ tGridIndexForDecomp ] = tGridIndices.size();
            tGridIndices.push_back( tGridIndexForDecomp );

            // add info about refinement
            moris_index tNumBoundaryRefs = tMaxTotalRefineLvlOnBspMeshes - tMaxUniformRefineLvlOnBspMeshes;
            tInitialRefinements.push_back( tMaxUniformRefineLvlOnBspMeshes );
            tBoundaryRefinements.push_back( tNumBoundaryRefs );

            // inform user that another grid is being used for decomposition
            MORIS_LOG( "Grid specified for decomposition is less refined than at least one B-spline mesh." );
            MORIS_LOG( "Creating grid #%i with %i uniform and %i interface refinements which will be used for decomposition.",
                    tGridIndexForDecomp,
                    tMaxUniformRefineLvlOnBspMeshes,
                    tNumBoundaryRefs );
        }

        // ------------------------------
        // Finalize the HMR parameter list

        // list the patterns with their initial refinements
        std::string sPatterns           = "";
        std::string sInitialRefinements = "";
        for ( uint iGrid = 0; iGrid < tGridIndices.size(); iGrid++ )
        {
            // add separators in string
            if ( iGrid > 0 )
            {
                sInitialRefinements += ",";
                sPatterns += ",";
            }

            // add it to the string
            sInitialRefinements += std::to_string( tInitialRefinements( iGrid ) );
            sPatterns += std::to_string( tGridIndices( iGrid ) );
        }

        // set initial refinement
        tHmrParamList.set( "initial_refinement_pattern", sPatterns );
        tHmrParamList.set( "initial_refinement", sInitialRefinements );

        // sort information about B-spline meshes
        std::string sLagrangeToBspline = "";
        std::string sBspPatterns       = "";
        std::string sBspPolyOrders     = "";
        for ( uint iBspMesh = 0; iBspMesh < tBspMeshIndices.size(); iBspMesh++ )
        {
            // get the number of refinements for this current grid
            moris_index tPattern   = tGridsForMeshes( iBspMesh );
            moris_index tPolyOrder = tPolyOrders( iBspMesh );

            // add it to the string
            if ( iBspMesh != 0 )
            {
                // add it to the string
                sLagrangeToBspline += ",";
                sBspPatterns += ",";
                sBspPolyOrders += ",";
            }

            sLagrangeToBspline += std::to_string( iBspMesh );
            sBspPatterns += std::to_string( tPattern );
            sBspPolyOrders += std::to_string( tPolyOrder );
        }

        // set the B-spline meshes in the parameter list
        tHmrParamList.set( "bspline_orders", sBspPolyOrders );
        tHmrParamList.set( "bspline_pattern", sBspPatterns );

        // create lagrange meshes for the decomposition grid ...
        std::string sLagrangeOrders   = std::to_string( tMaxPolyOrder );
        std::string sLagrangePatterns = std::to_string( tGridIndexForDecomp );
        mGenNumRefinements            = { (uint)tBoundaryRefinements( tGridIndexForDecomp ) };

        // create dummy lagrange meshes using the other grids to trigger geometric refinement through these
        std::string sLagrangeToBsplineAddOn = "";
        mGenRefineMeshIndices               = { 0 };
        uint tLagMeshCounter                = 1;    // start at 1 as Lagrange mesh index 0 is already occupied by the decomposition grid
        for ( uint iGrid = 0; iGrid < tGridIndices.size(); iGrid++ )
        {
            moris_index tPatternIndex       = tGridIndices( iGrid );
            moris_index tNumBdRefsOnPattern = tBoundaryRefinements( iGrid );

            // don't create a dummy lagrange mesh for the grid specified for decomposition
            if ( tPatternIndex != tGridIndexForDecomp )
            {
                sLagrangeOrders += "," + std::to_string( tMaxPolyOrder );
                sLagrangePatterns += ( "," + std::to_string( tPatternIndex ) );
                sLagrangeToBsplineAddOn += ";-1";

                // store boundary refinements associated with this pattern to later feed through the GEN parameter list
                mGenRefineMeshIndices.push_back( tLagMeshCounter );
                mGenNumRefinements.push_back( tNumBdRefsOnPattern );
                tLagMeshCounter++;
            }
        }

        // all parameters that have been retrieved up to this point (set pattern information and Lagrange mesh)
        tHmrParamList.set( "lagrange_orders", sLagrangeOrders );
        tHmrParamList.set( "lagrange_pattern", sLagrangePatterns );
        tHmrParamList.set( "lagrange_to_bspline", sLagrangeToBspline + sLagrangeToBsplineAddOn );

        // set the refinement buffer such that the refinement actually has an effect on the highest order B-spline mesh
        int tRefinementBuffer = std::max( (int)tMaxPolyOrder - 1, 1 );
        tHmrParamList.set( "refinement_buffer", tRefinementBuffer );

    }    // end function: Library_IO_Meshgen::load_HMR_parameters_from_xml()

    // -----------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_XTK_parameters_from_xml(
            std::string const & aXtkPath,
            std::string const & aHmrPath )
    {
        // quick access to the parameter list
        Parameter_List& tXtkParamList = mParameterLists( (uint)( Module_Type::XTK ) )( 0 )( 0 );
        Parameter_List& tHmrParamList = mParameterLists( (uint)( Module_Type::HMR ) )( 0 )( 0 );

        // turn on SPG based enrichment to make sure
        tXtkParamList.set( "use_SPG_based_enrichment", true );

        // enriched mesh indices
        std::string sLagrangeToBspline  = tHmrParamList.get< std::string >( "lagrange_to_bspline" );
        std::string sBsplineMeshIndices = sLagrangeToBspline.substr( 0, sLagrangeToBspline.find( ';' ) );
        tXtkParamList.set( "enrich_mesh_indices", sBsplineMeshIndices );

        // get whether to triangulate all
        bool tTriangulateAll = false;
        mXmlReader->get( aXtkPath + ".TriangulateAllFgElems", tTriangulateAll, false );
        tXtkParamList.set( "triangulate_all", tTriangulateAll );

        // check that boundary refinement is not requested when triangulating all
        if ( tTriangulateAll )
        {
            // get the number of geometric refinements around geometric boundaries
            moris_index tNumBoundaryRefinements = 0;
            mXmlReader->get( aXtkPath + ".InterfaceRefinements", tNumBoundaryRefinements, 0 );

            MORIS_ERROR(
                    tNumBoundaryRefinements == 0,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Triangulation of all elements and boundary refinement at the same time are not supported yet due to hanging nodes." );
        }

        // option to output Lagrange meshes
        bool tOutputDecompGrid = false;
        mXmlReader->get( aXtkPath + ".OutputDecompositionGrid", tOutputDecompGrid, false );
        if ( tOutputDecompGrid )
        {
            tHmrParamList.set( "write_lagrange_output_mesh_to_exodus", "Decomposition_Grid.exo" );
        }

        // get the number of refinements of the Lagrange mesh at the boundary
        bool tUseBasisExtensions = false;
        mXmlReader->get( aXtkPath + ".UseCutBasisAgglomeration", tUseBasisExtensions, false );
        tXtkParamList.set( "activate_basis_agglomeration", tUseBasisExtensions );

        // get foreground mesh order
        uint tFgElemPolyOrder = 1;
        mXmlReader->get( aXtkPath + ".FgPolynomialOrder", tFgElemPolyOrder, uint( 1 ) );
        MORIS_ERROR( tFgElemPolyOrder == 1 || tFgElemPolyOrder == 2,
                "Library_IO_Meshgen::load_parameters_from_xml() - Currently only supporting foreground polynomial orders 1 and 2." );
        tXtkParamList.set( "ig_element_order", tFgElemPolyOrder );

        // check whether T-matrix output has been requested/suppressed
        bool tOutputTmats = "";
        mXmlReader->get( aXtkPath + ".OutputExtractionOperators", tOutputTmats, true );
        tXtkParamList.set( "only_generate_xtk_temp", !tOutputTmats );

        // check which T-matrix outputs have been requested
        std::string tTmatOutputFormats = "";
        mXmlReader->get( aXtkPath + ".ExtractionOperatorFormat", tTmatOutputFormats, std::string( "" ) );
        bool tOutputElemental = ( tTmatOutputFormats.find( "Elemental" ) != std::string::npos );
        bool tOutputGlobal    = ( tTmatOutputFormats.find( "Global" ) != std::string::npos );

        // if non have been specified just output elemental
        if ( !tOutputElemental && !tOutputElemental )
        {
            tOutputElemental = true;
        }

        // set the T-matrix outputs
        if ( tOutputGlobal )
        {
            tXtkParamList.set( "global_T_matrix_output_file", "Global_Extraction_Operators" );
        }
        if ( tOutputElemental )
        {
            tXtkParamList.set( "elemental_T_matrix_output_file", "Elemental_Extraction_Operators" );
        }

    }    // end function: Library_IO_Meshgen::load_XTK_parameters_from_xml()

    // -----------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_GEN_parameters_from_xml(
            std::string const & aGenPath,
            std::string const & aHmrPath,
            std::string const & aXtkPath )
    {
        // quick access to the parameter list
        Module_Parameter_Lists& tGenParamList = mParameterLists( (uint)( Module_Type::GEN ) );

        // path to the meshes
        std::string tGeometryNodeName = "Geometry";

        // see how many grids are specified in the input file
        uint tNumGeometries = mXmlReader->count_keys_in_subtree( aGenPath, tGeometryNodeName );

        // get the intersection mode
        bool tUseMultiLinearIntersections = true;
        mXmlReader->get( aGenPath + ".UseMultiLinearIntersections", tUseMultiLinearIntersections, true );

        // get the number of geometric refinements around geometric boundaries
        moris_index tNumBoundaryRefinements = 0;
        mXmlReader->get( aXtkPath + ".InterfaceRefinements", tNumBoundaryRefinements, 0 );

        // go over geometries and load their respective parameters
        for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
        {
            // move current geometry into the buffer
            mXmlReader->copy_subtree_into_buffer( aGenPath, tGeometryNodeName, iGeom );

            // get the type to geometry defined
            std::string tGeomType = mXmlReader->get_attribute_from_buffer( "field_type", std::string( "" ) );

            // check that a type is defined
            MORIS_ERROR( tGeomType != "",
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "No type defined for this geometry. For all geometries an attribute 'type' needs to be defined." );

            // create parameter lists according to geometry type
            if ( tGeomType == "pre_defined" )
            {
                // find the 'geom' attribute
                std::string tPreDefGeom = mXmlReader->get_attribute_from_buffer( "geom", std::string( "" ) );
                MORIS_ERROR( tPreDefGeom != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "All pre-defined geometries must have an attribute 'geom' specified. Supported Options are 'plane' and 'circle'." );

                // -------------------------------- //
                // PLANE

                if ( tPreDefGeom == "line" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::LINE );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All pre-defined geometries must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == mNumSpatialDims || tPointMat.n_cols() == mNumSpatialDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector does not match number of spatial dimensions" );

                    // get the normal
                    std::string tNormal = "";
                    mXmlReader->get_from_buffer( "Normal", tNormal, std::string( "" ) );
                    MORIS_ERROR( tNormal != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All planes must have a parameter 'Normal' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tNormalMat;
                    moris::string_to_matrix( tNormal, tNormalMat );
                    MORIS_ERROR( tNormalMat.numel() == mNumSpatialDims || tNormalMat.n_cols() == mNumSpatialDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Normal' vector does not match number of spatial dimensions for the 'plane'." );

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector  = string_to_vector< real >( tPoint );
                    Vector< real > tNormalVector = string_to_vector< real >( tNormal );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "normal_x", tNormalVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "normal_y", tNormalVector( 1 ) );
                }

                // -------------------------------- //
                // CIRCLE

                else if ( tPreDefGeom == "circle" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::CIRCLE );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'circle' must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == 2 || tPointMat.n_cols() == 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector for a 'circle' needs to be two." );

                    // get the radius
                    std::string tRadius = "";
                    mXmlReader->get_from_buffer( "Radius", tRadius, std::string( "" ) );
                    MORIS_ERROR( tRadius != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All planes must have a parameter 'Radius' specified of format e.g.: '5.6'" );

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector = string_to_vector< real >( tPoint );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "radius", std::stod( tRadius ) );
                }

                // -------------------------------- //
                // ELLIPSE

                else if ( tPreDefGeom == "ellipse" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::SUPERELLIPSE );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipse' must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == 2 || tPointMat.n_cols() == 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector for a 'ellipse' needs to be two." );

                    // get the half diameters
                    std::string tSemiDiameters = "";
                    mXmlReader->get_from_buffer( "SemiDiameters", tSemiDiameters, std::string( "" ) );
                    MORIS_ERROR( tSemiDiameters != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipse' must have a parameter 'SemiDiameters' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tSemiDiametersMat;
                    moris::string_to_matrix( tSemiDiameters, tSemiDiametersMat );
                    MORIS_ERROR( tSemiDiametersMat.numel() == 2 || tSemiDiametersMat.n_cols() == 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'SemiDiameters' vector for a 'ellipse' needs to be two." );

                    // get the exponent, default to 2.0 if not specified
                    std::string tExponent = "2.0";
                    mXmlReader->get_from_buffer( "Exponent", tExponent, std::string( "2.0" ) );
                    if ( tExponent == "" ) { tExponent = "2.0"; };

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector        = string_to_vector< real >( tPoint );
                    Vector< real > tSemidiameterVector = string_to_vector< real >( tSemiDiameters );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "semidiameter_x", tSemidiameterVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "semidiameter_y", tSemidiameterVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "exponent", std::stod( tExponent ) );
                }

                else if ( tPreDefGeom == "plane" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::PLANE );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All pre-defined geometries must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == mNumSpatialDims || tPointMat.n_cols() == mNumSpatialDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector does not match number of spatial dimensions" );

                    // get the normal
                    std::string tNormal = "";
                    mXmlReader->get_from_buffer( "Normal", tNormal, std::string( "" ) );
                    MORIS_ERROR( tNormal != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All planes must have a parameter 'Normal' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tNormalMat;
                    moris::string_to_matrix( tNormal, tNormalMat );
                    MORIS_ERROR( tNormalMat.numel() == mNumSpatialDims || tNormalMat.n_cols() == mNumSpatialDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Normal' vector does not match number of spatial dimensions for the 'plane'." );

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector  = string_to_vector< real >( tPoint );
                    Vector< real > tNormalVector = string_to_vector< real >( tNormal );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_z", tPointVector( 2 ) );
                    tGenParamList( 1 )( iGeom ).set( "normal_x", tNormalVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "normal_y", tNormalVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "normal_z", tNormalVector( 2 ) );
                }

                // -------------------------------- //
                // SPHERE

                else if ( tPreDefGeom == "sphere" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::SPHERE );

                    // check dimensionality
                    MORIS_ERROR( mNumSpatialDims != 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "A 'sphere' can only be defined in 3D, but mesh is 2D." );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'circle' must have a parameter 'Point' specified of format e.g.: '1.2,3.4,5.6'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == 3 || tPointMat.n_cols() == 3,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector for a 'sphere' needs to be three." );

                    // get the radius
                    std::string tRadius = "";
                    mXmlReader->get_from_buffer( "Radius", tRadius, std::string( "" ) );
                    MORIS_ERROR( tRadius != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All planes must have a parameter 'Radius' specified of format e.g.: '5.6'" );

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector = string_to_vector< real >( tPoint );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_z", tPointVector( 2 ) );
                    tGenParamList( 1 )( iGeom ).set( "radius", std::stod( tRadius ) );
                }

                // -------------------------------- //
                // ELLIPSE

                else if ( tPreDefGeom == "ellipsoid" )
                {
                    // Create geometry parameter list
                    tGenParamList( 1 ).add_parameter_list( gen::Field_Type::SUPERELLIPSOID );

                    // check dimensionality
                    MORIS_ERROR( mNumSpatialDims != 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "A 'sphere' can only be defined in 3D, but mesh is 2D." );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipsoid' must have a parameter 'Point' specified of format e.g.: '1.2,3.4,5.6'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_matrix( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == 3 || tPointMat.n_cols() == 3,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector for a 'ellipsoid' needs to be three." );

                    // get the half diameters
                    std::string tSemiDiameters = "";
                    mXmlReader->get_from_buffer( "Point", tSemiDiameters, std::string( "" ) );
                    MORIS_ERROR( tSemiDiameters != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipsoid' must have a parameter 'SemiDiameters' specified of format e.g.: '1.2,3.4,5.6'" );
                    Matrix< DDRMat > tSemiDiametersMat;
                    moris::string_to_matrix( tSemiDiameters, tSemiDiametersMat );
                    MORIS_ERROR( tSemiDiametersMat.numel() == 3 || tSemiDiametersMat.n_cols() == 3,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'SemiDiameters' vector for a 'ellipsoid' needs to be three." );

                    // get the exponent, default to 2.0 if not specified
                    std::string tExponent = "2.0";
                    mXmlReader->get_from_buffer( "Exponent", tExponent, std::string( "2.0" ) );
                    if ( tExponent == "" ) { tExponent = "2.0"; };

                    // set the parameters in the GEN parameter list
                    Vector< real > tPointVector        = string_to_vector< real >( tPoint );
                    Vector< real > tSemidiameterVector = string_to_vector< real >( tSemiDiameters );
                    tGenParamList( 1 )( iGeom ).set( "center_x", tPointVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_y", tPointVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "center_z", tPointVector( 2 ) );
                    tGenParamList( 1 )( iGeom ).set( "semidiameter_x", tSemidiameterVector( 0 ) );
                    tGenParamList( 1 )( iGeom ).set( "semidiameter_y", tSemidiameterVector( 1 ) );
                    tGenParamList( 1 )( iGeom ).set( "semidiameter_z", tSemidiameterVector( 2 ) );
                    tGenParamList( 1 )( iGeom ).set( "exponent", std::stod( tExponent ) );
                }

                // -------------------------------- //
                // UNKNOWN GEOMETRY

                else
                {
                    MORIS_ERROR( false,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Unknown pre-defined geometry. Supported Options are 'plane', 'circle', 'sphere', 'ellipse', and 'ellipsoid'." );
                }

                // set parameters that are independent of the particular geometry used
                tGenParamList( 1 )( iGeom ).set( "number_of_refinements", mGenNumRefinements );
                tGenParamList( 1 )( iGeom ).set( "refinement_mesh_index", mGenRefineMeshIndices );
                tGenParamList( 1 )( iGeom ).set( "use_multilinear_interpolation", tUseMultiLinearIntersections );
                tGenParamList( 1 )( iGeom ).set( "discretization_mesh_index", -1 );

            }    // end if: pre-defined geometry

            else if ( tGeomType == "image_file" )
            {
                // initialize with the image sdf default parameter list
                tGenParamList( 1 ).add_parameter_list( gen::Field_Type::SIGNED_DISTANCE_IMAGE );
                tGenParamList( 1 )( iGeom ).set( "number_of_refinements", mGenNumRefinements );
                tGenParamList( 1 )( iGeom ).set( "refinement_mesh_index", mGenRefineMeshIndices );
                tGenParamList( 1 )( iGeom ).set( "use_multilinear_interpolation", tUseMultiLinearIntersections );
                tGenParamList( 1 )( iGeom ).set( "discretization_mesh_index", -1 );

                // get the file name and check it for validity
                std::string tFileName = "";
                mXmlReader->get_from_buffer( "FileName", tFileName, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No image file name provided for geometry #%i. Please provide an .hdf5 image file using the tag 'FileName' ",
                        iGeom );
                MORIS_ERROR( this->string_ends_with( tFileName, ".h5" ) || this->string_ends_with( tFileName, ".hdf5" ),
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "File provided is not an .h5 or .hdf5 file. Image data needs to be provided in this format." );
                tGenParamList( 1 )( iGeom ).set( "image_file", tFileName );

                // get the image dimensions
                std::string tImageDimStr = "";
                mXmlReader->get_from_buffer( "ImageDimensions", tImageDimStr, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No image dimensions provided for geometry #%i. "
                        "Please provide image dimensions using the tag 'ImageDimensions' in the the format e.g. '1.2,3.3' ",
                        iGeom );
                Vector< real > tImageDimVec;
                moris::string_to_vector( tImageDimStr, tImageDimVec );
                MORIS_ERROR( tImageDimVec.size() == mNumSpatialDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ImageDimensions' vector does not match number of spatial dimensions" );
                tGenParamList( 1 )( iGeom ).set( "image_dimensions", tImageDimVec );

                // get the image offset
                std::string tImageOffsetStr = "";
                mXmlReader->get_from_buffer( "ImageOrigin", tImageOffsetStr, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No image origin provided for geometry #%i. "
                        "Please provide image origin/offset using the tag 'ImageOrigin' in the the format e.g. '1.2,3.3' ",
                        iGeom );
                Vector< real > tImageOffsetVec;
                moris::string_to_vector( tImageOffsetStr, tImageOffsetVec );
                MORIS_ERROR( tImageOffsetVec.size() == mNumSpatialDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ImageOrigin' vector does not match number of spatial dimensions" );
                tGenParamList( 1 )( iGeom ).set( "image_offset", tImageOffsetVec );

            }    // end if: geometry from image file

            else if ( tGeomType == "object_file" )
            {
                // initialize with the sdf field default parameter list
                tGenParamList( 1 ).add_parameter_list( gen::Field_Type::SIGNED_DISTANCE_OBJECT );
                tGenParamList( 1 )( iGeom ).set( "use_multilinear_interpolation", false );
                // tGenParamList( 1 )( iGeom ).set( "discretization_mesh_index", -1 );

                // FIXME: this functionality needs to get added
                // tGenParamList( 1 )( iGeom ).set( "number_of_refinements", mGenNumRefinements );
                // tGenParamList( 1 )( iGeom ).set( "refinement_mesh_index", mGenRefineMeshIndices );

                // get the file name and check it for validity
                std::string tFileName = "";
                mXmlReader->get_from_buffer( "FileName", tFileName, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No object file name provided for geometry #%i. Please provide an .obj file using the tag 'FileName' ",
                        iGeom );
                MORIS_ERROR( this->string_ends_with( tFileName, ".obj" ),
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "File provided is not an .obj file. STL data needs to be provided in this format." );
                tGenParamList( 1 )( iGeom ).set( "sdf_object_path", tFileName );

                // get the object offset
                std::string tObjectOffsetStr = "";
                mXmlReader->get_from_buffer( "ObjectOrigin", tObjectOffsetStr, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No object origin provided for geometry #%i. "
                        "Please provide an origin/offset using the tag 'ObjectOrigin' in the the format e.g. '1.2,3.3' ",
                        iGeom );
                Vector< real > tObjectOffsetVec;
                moris::string_to_vector( tObjectOffsetStr, tObjectOffsetVec );
                MORIS_ERROR( tObjectOffsetVec.size() == mNumSpatialDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ObjectOrigin' vector for geometry %i does not match number of spatial dimensions.",
                        iGeom );
                tGenParamList( 1 )( iGeom ).set( "sdf_object_offset", tObjectOffsetVec );

                // get an offset in the sign distance value if input is provided
                double tSdfShift = 0.0;
                mXmlReader->get_from_buffer( "SdfShift", tSdfShift, 0.0 );
                tGenParamList( 1 )( iGeom ).set( "sdf_shift", tSdfShift );

            }    // end if: geometry from image file
            else
            {
                MORIS_ERROR( false,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Geometry type '%s' unknown. Currently supported options are: 'pre_defined', 'image_file' and 'object_file'.",
                        tGeomType.c_str() );
            }
        }    // end for: loop over geometries in input file

        // ------------------------------
        // Phase map

        // get the string from the input file
        std::string tPhaseMapString = "";
        mXmlReader->get( aGenPath + ".PhaseMap", tPhaseMapString, std::string( "" ) );
        bool tPhaseMapSpecified = ( tPhaseMapString != "" );

        // if a phase map has been specified, convert it to a matrix
        Matrix< DDUMat > tPhaseMap;
        moris::string_to_matrix( tPhaseMapString, tPhaseMap );

        if ( tPhaseMapSpecified )
        {
            // make sure the input is set correctly
            MORIS_ERROR( tPhaseMap.n_cols() == 2,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Format of the phase map incorrect. For correct use always map one phase index to one material index, e.g.: 0,1;1,1;2,0" );

            // initialize phase table
            uint             tNumPhases = std::pow( 2, tNumGeometries );
            Matrix< DDUMat > tPhaseTable( tNumPhases, 1, 0 );

            // fill phase table with default values
            for ( uint i = 0; i < tNumPhases; i++ )
            {
                tPhaseTable( i ) = i;
            }

            // use the map to correct the table
            for ( uint i = 0; i < tPhaseMap.n_rows(); i++ )
            {
                uint tPhaseToMapFrom           = tPhaseMap( i, 0 );
                uint tMaterialToMapTo          = tPhaseMap( i, 1 );
                tPhaseTable( tPhaseToMapFrom ) = tMaterialToMapTo;
            }

            // convert the phase table back to a string
            std::string tPhaseTableString = ios::stringify( tPhaseTable );

            // set the phase table in the parameter list
            tGenParamList( 0 )( 0 ).set( "phase_table", tPhaseTableString );
        }

    }    // end function: Library_IO_Meshgen::load_GEN_parameters_from_xml()

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_all_standard_parameters()
    {
        // go over all modules and check that
        for ( uint iModule = 0; iModule < (uint)( Module_Type::END_ENUM ); iModule++ )
        {
            // get the current module
            Module_Type tParamListType = (Module_Type)( iModule );

            // fill the parameter list entry with the standard parameters
            this->create_standard_parameter_list_for_module( tParamListType, mParameterLists( iModule ) );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_parameter_list_for_module(
            Module_Type             aParamListType,
            Module_Parameter_Lists& aParameterList )
    {
        switch ( aParamListType )
        {
            case Module_Type::OPT:
                this->create_standard_OPT_parameter_list( aParameterList );
                break;

            case Module_Type::HMR:
                this->create_standard_HMR_parameter_list( aParameterList );
                break;

            case Module_Type::XTK:
                this->create_standard_XTK_parameter_list( aParameterList );
                break;

            case Module_Type::GEN:
                this->create_standard_GEN_parameter_list( aParameterList );
                break;

            case Module_Type::END_ENUM:
                MORIS_ERROR( false,
                        "Library_IO_Meshgen::create_standard_parameter_list_for_module() - "
                        "No standard library defined for module UNDEFINED" );
                break;

            // create an empty parameter list for modules that are not needed
            default:
                aParameterList = Module_Parameter_Lists( Module_Type::END_ENUM );
                break;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_OPT_parameter_list( Module_Parameter_Lists& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_XTK_parameter_list( Module_Parameter_Lists& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );

        // enrichment
        aParameterList( 0 )( 0 ).set( "enrich", true );
        aParameterList( 0 )( 0 ).set( "basis_rank", "bspline" );

        // output
        aParameterList( 0 )( 0 ).set( "output_file", "foreground_mesh.exo" );
        aParameterList( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterList( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterList( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_HMR_parameter_list( Module_Parameter_Lists& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        // Lagrange mesh is always 0
        aParameterList( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_GEN_parameter_list( Module_Parameter_Lists& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
