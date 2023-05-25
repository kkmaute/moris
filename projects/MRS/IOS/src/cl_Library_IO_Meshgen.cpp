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
        // set the type of this library
        mLibraryType = Library_Type::MESHGEN;

        // list of supported parameter list types
        mSupportedParamListTypes = { Parameter_List_Type::HMR, Parameter_List_Type::XTK, Parameter_List_Type::GEN };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Meshgen::~Library_IO_Meshgen()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::finalize()
    {
        // check that an .xml input file has been specified
        MORIS_ERROR( mXmlParserIsInitialized,
                "Library_IO_Meshgen::finalize() - No .xml input file has been specified. "
                "This is required for the mesh generation workflow." );

        // load the standard parameters into the member variables
        this->load_all_standard_parameters();

        // if an .so file has been parsed, first use its parameters (if any were defined in it) to overwrite or add to the standard parameters
        if ( mSoLibIsInitialized )
        {
            this->load_parameters_from_shared_object_library();
        }

        // load parameters from xml, overwrites parameters specified in either the standard parameters or an .so file if parsed
        if ( mXmlParserIsInitialized )
        {
            this->load_parameters_from_xml();
        }

        // check the parameters for validity
        // TODO: this->check_parameters();

        // mark this library as finalized and lock it from modification
        mLibraryIsFinalized = true;

        // print receipt of the finalized library
        this->print_parameter_receipt( "./Parameter_Receipt.xml" );    // TODO: the file name and location should be user define-able
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

        // quick access to the various parameter lists
        ParameterList&       tHmrParamList = mParameterLists( (uint)( Parameter_List_Type::HMR ) )( 0 )( 0 );
        ParameterList&       tXtkParamList = mParameterLists( (uint)( Parameter_List_Type::XTK ) )( 0 )( 0 );
        ModuleParameterList& tGenParamList = mParameterLists( (uint)( Parameter_List_Type::GEN ) );

        // ------------------------------
        // Base grid

        // get the base grid parameters
        std::string tBaseGridPath = tHmrPath + ".BaseGrid";
        std::string tBaseGridSize;
        std::string tDomainDimensions;
        std::string tBaseGridOrigin;
        mXmlReader->get( tBaseGridPath + ".Size", tBaseGridSize, std::string( "" ) );
        mXmlReader->get( tBaseGridPath + ".Dimensions", tDomainDimensions, std::string( "" ) );
        mXmlReader->get( tBaseGridPath + ".Origin", tBaseGridOrigin, std::string( "" ) );

        // check the number of spatial dimensions
        Matrix< DDUMat > tBaseGridMat;
        moris::string_to_mat( tBaseGridSize, tBaseGridMat );
        uint tNumDims = tBaseGridMat.numel();

        Matrix< DDRMat > tDomainDimsMat;
        moris::string_to_mat( tDomainDimensions, tDomainDimsMat );
        uint tNumDimsDimensions = tDomainDimsMat.numel();

        Matrix< DDRMat > tDomainOriginMat;
        moris::string_to_mat( tBaseGridOrigin, tDomainOriginMat );
        uint tNumDimsOrigin = tDomainOriginMat.numel();

        // perform some checks on the user inputs
        MORIS_ERROR( tNumDims == 2 || tNumDims == 3,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid specified must have 2 or 3 dimensions." );
        MORIS_ERROR( tNumDims == tNumDimsDimensions,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid 'Dimensions' inconsistent with 'Size'." );
        MORIS_ERROR( tNumDims == tNumDimsOrigin,
                "Library_IO_Meshgen::load_parameters_from_xml() - Base grid 'Origin' inconsistent with 'Size'." );

        // generate the correct side set names for this number of dimensions
        std::string tDomainSideSets = ( tNumDims == 3 ) ? "1,2,3,4,5,6" : "1,2,3,4";

        // set the parameters in the parameter lists
        tHmrParamList.set( "number_of_elements_per_dimension", tBaseGridSize );
        tHmrParamList.set( "domain_dimensions", tDomainDimensions );
        tHmrParamList.set( "domain_offset", tBaseGridOrigin );
        tHmrParamList.set( "domain_sidesets", tDomainSideSets );

        // ------------------------------
        // grids and refinements

        // path to the grids
        std::string tGridsPath    = tHmrPath + ".MeshGrids";
        std::string tGridNodeName = "MeshGrid";

        // see how many grids are specified in the input file
        uint tNumGridsSpecified = mXmlReader->count_keys_in_subtree( tGridsPath, tGridNodeName );

        // initialize a map that stores in what order the grids are defined
        std::unordered_map< moris_index, uint > tGridIndexMap;
        Cell< moris_index >                     tRefinements( tNumGridsSpecified, 0 );
        moris_index                             tMaxGridIndex = -1;

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

            // get and store the number of refinements
            moris_index& tNumRefines = tRefinements( iMeshGrid );
            mXmlReader->get_from_buffer( "InitialRefinements", tNumRefines, moris_index( 0 ) );
        }

        // make sure there aren't any mesh indices left unidentified
        MORIS_ERROR( tMaxGridIndex == (moris_index)(tNumGridsSpecified)-1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Some grid indices have not been defined. The maximum grid index used is %i, but only %i grids are defined.",
                tMaxGridIndex,
                tNumGridsSpecified );

        // ------------------------------
        // B-spline meshes

        // path to the meshes
        std::string tBspMeshesPath   = tHmrPath + ".BsplineMeshes";
        std::string tBspMeshNodeName = "BsplineMesh";

        // see how many grids are specified in the input file
        uint tNumBspMeshesSpecified = mXmlReader->count_keys_in_subtree( tBspMeshesPath, tBspMeshNodeName );

        // initialize a map that stores in what order the grids are defined
        std::unordered_map< moris_index, uint > tBspMeshIndexMap;
        Cell< moris_index >                     tGridsForMeshes( tNumBspMeshesSpecified, -1 );
        Cell< moris_index >                     tPolyOrders( tNumBspMeshesSpecified, -1 );
        moris_index                             tMaxBspMeshIndex = -1;
        moris_index                             tMaxPolyOrder    = -1;

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
                    "B-spline mesh index %i defined twice. Please use unique indices.",
                    tBspMeshIndex );

            // add to map
            tBspMeshIndexMap[ tBspMeshIndex ] = iBspMesh;

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
        moris_index tGridForDecomp = 0;
        mXmlReader->get( tXtkPath + ".DecompositionGrid", tGridForDecomp, 0 );
        MORIS_ERROR( tGridForDecomp > -1 && tGridForDecomp < tMaxGridIndex + 1,
                "Library_IO_Meshgen::load_parameters_from_xml() - "
                "Trying to use mesh grid index %i for decomposition. But only grid indices 0 to %i are defined.",
                tGridForDecomp,
                tMaxGridIndex );

        // get the number of refinements specified for the decomp mesh
        auto tDecompGridIter = tGridIndexMap.find( tGridForDecomp );
        MORIS_ERROR( tDecompGridIter != tGridIndexMap.end(),
                "Library_IO_Meshgen::load_parameters_from_xml() - Grid index for decomposition not found in map. Something went wrong." );
        moris_index tNumInitRefinesDecompGrid = tRefinements( tDecompGridIter->second );

        // get the number of refinements of the Lagrange mesh at the boundary
        moris_index tNumBoundaryRefinements = 0;
        mXmlReader->get( tXtkPath + ".InterfaceRefinements", tNumBoundaryRefinements, 0 );

        // create additional pattern for the lagrange mesh
        moris_index tLagrangePattern = tMaxGridIndex + 1;

        // ------------------------------
        // Finalize the HMR parameter list

        // list the patterns with their initial refinements
        std::string tPatterns           = "";
        std::string tInitialRefinements = "";
        for ( uint iGrid = 0; iGrid < (uint)tMaxGridIndex + 1; iGrid++ )
        {
            // get the postion in the lists
            auto tIter = tGridIndexMap.find( iGrid );
            MORIS_ERROR( tIter != tGridIndexMap.end(), "Library_IO_Meshgen::load_parameters_from_xml() - Grid index not found in map. Something went wrong." );
            uint tPos = tIter->second;

            // get the number of refinements for this current grid
            moris_index tNumRefines = tRefinements( tPos );

            // add separators in string
            if ( iGrid > 0 )
            {
                tInitialRefinements += ",";
                tPatterns += ",";
            }

            // add it to the string
            tInitialRefinements += std::to_string( tNumRefines );
            tPatterns += std::to_string( iGrid );
        }

        // append the pattern which will be used for the lagrange mesh
        tInitialRefinements += "," + std::to_string( tNumInitRefinesDecompGrid );
        tPatterns += "," + std::to_string( tMaxGridIndex + 1 );

        // all parameters that have been retrieved up to this point (set pattern information and Lagrange mesh)
        tHmrParamList.set( "initial_refinement_pattern", tPatterns );
        tHmrParamList.set( "initial_refinement", tInitialRefinements );
        tHmrParamList.set( "lagrange_orders", std::to_string( tMaxPolyOrder ) );
        tHmrParamList.set( "lagrange_pattern", std::to_string( tLagrangePattern ) );

        // sort information about B-spline meshes
        std::string tLagrangeToBspline = "";
        std::string tBspPatterns       = "";
        std::string tBspPolyOrders     = "";
        for ( uint iBspMesh = 0; iBspMesh < (uint)tMaxBspMeshIndex + 1; iBspMesh++ )
        {
            // get the postion in the lists
            auto tIter = tBspMeshIndexMap.find( iBspMesh );
            MORIS_ERROR( tIter != tGridIndexMap.end(),
                    "Library_IO_Meshgen::load_parameters_from_xml() - B-spline mesh index not found in map. Something went wrong." );
            uint tPos = tIter->second;

            // get the number of refinements for this current grid
            moris_index tPattern   = tGridsForMeshes( tPos );
            moris_index tPolyOrder = tPolyOrders( tPos );

            // add it to the string
            if ( iBspMesh != 0 )
            {
                // add it to the string
                tLagrangeToBspline += ",";
                tBspPatterns += ",";
                tBspPolyOrders += ",";
            }
            tLagrangeToBspline += std::to_string( iBspMesh );
            tBspPatterns += std::to_string( tPattern );
            tBspPolyOrders += std::to_string( tPolyOrder );
        }

        // set the B-spline meshes in the parameter list
        tHmrParamList.set( "bspline_orders", tBspPolyOrders );
        tHmrParamList.set( "bspline_pattern", tBspPatterns );
        tHmrParamList.set( "lagrange_to_bspline", tLagrangeToBspline );

        // ------------------------------
        // XTK parameters

        // turn on SPG based enrichment to make sure 
        tXtkParamList.set( "use_SPG_based_enrichment", true );

        // enriched mesh indices
        tXtkParamList.set( "enrich_mesh_indices", tLagrangeToBspline );

        // get whether to triangulate all
        bool tTriangulateAll = false;
        mXmlReader->get( tXtkPath + ".TriangulateAllFgElems", tTriangulateAll, bool( false ) );
        tXtkParamList.set( "triangulate_all", tTriangulateAll );

        // check that boundary refinement is not requested when triangulating all
        if ( tTriangulateAll )
        {
            MORIS_ERROR( tNumBoundaryRefinements == 0,
                    "Library_IO_Meshgen::load_parameters_from_xml() - "
                    "Triangulation of all elements and boundary refinement at the same time are not supported yet due to hanging nodes." );
        }

        // get foreground mesh order
        uint tFgElemPolyOrder = 1;
        mXmlReader->get( tXtkPath + ".FgPolynomialOrder", tFgElemPolyOrder, uint( 1 ) );
        MORIS_ERROR( tFgElemPolyOrder == 1 || tFgElemPolyOrder == 2,
                "Library_IO_Meshgen::load_parameters_from_xml() - Currently only supporting foreground polynomial orders 1 and 2." );
        tXtkParamList.set( "ig_element_order", tFgElemPolyOrder );

        // check whether T-matrix output has been requested/suppressed 
        bool tOutputTmats = "";
        mXmlReader->get( tXtkPath + ".OutputExtractionOperators", tOutputTmats, bool( true ) );
        tXtkParamList.set( "only_generate_xtk_temp", !tOutputTmats );

        // check which T-matrix outputs have been requested
        std::string tTmatOutputFormats = "";
        mXmlReader->get( tXtkPath + ".ExtractionOperatorFormat", tTmatOutputFormats, std::string( "" ) );
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

        // ------------------------------
        // Geometries

        // path to the meshes
        std::string tGeometryNodeName = "Geometry";

        // see how many grids are specified in the input file
        uint tNumGeometries = mXmlReader->count_keys_in_subtree( tGenPath, tGeometryNodeName );

        // resize the parameter list correctly
        tGenParamList( 1 ).resize( tNumGeometries );

        // go over geometries and load their respective parameters
        for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
        {
            // load the standard geometry parameter regardless of type
            tGenParamList( 1 )( iGeom ) = prm::create_geometry_parameter_list();
            tGenParamList( 1 )( iGeom ).set( "number_of_refinements", std::to_string( tNumBoundaryRefinements ) );
            tGenParamList( 1 )( iGeom ).set( "refinement_mesh_index", "0" );

            // move current geometry into the buffer
            mXmlReader->copy_subtree_into_buffer( tGenPath, tGeometryNodeName, iGeom );

            // get the type to geometry defined
            std::string tGeomType = mXmlReader->get_attribute_from_buffer( "type", std::string( "" ) );

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

                if ( tPreDefGeom == "plane" )
                {
                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All pre-defined geometries must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_mat( tPoint, tPointMat );
                    MORIS_ERROR( tPointMat.numel() == tNumDims || tPointMat.n_cols() == tNumDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Point' vector does not match number of spatial dimensions" );

                    // get the normal
                    std::string tNormal = "";
                    mXmlReader->get_from_buffer( "Normal", tNormal, std::string( "" ) );
                    MORIS_ERROR( tNormal != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "All planes must have a parameter 'Normal' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tNormalMat;
                    moris::string_to_mat( tNormal, tNormalMat );
                    MORIS_ERROR( tNormalMat.numel() == tNumDims || tNormalMat.n_cols() == tNumDims,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'Normal' vector does not match number of spatial dimensions for the 'plane'." );

                    // set the parameters in the GEN parameter list
                    tGenParamList( 1 )( iGeom ).set( "type", "plane" );
                    tGenParamList( 1 )( iGeom ).set( "constant_parameters", tPoint + "," + tNormal );
                }

                // -------------------------------- //
                // CIRCLE

                else if ( tPreDefGeom == "circle" )
                {
                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'circle' must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_mat( tPoint, tPointMat );
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
                    tGenParamList( 1 )( iGeom ).set( "type", "circle" );
                    tGenParamList( 1 )( iGeom ).set( "constant_parameters", tPoint + "," + tRadius );
                }

                // -------------------------------- //
                // ELLIPSE

                else if ( tPreDefGeom == "ellipse" )
                {
                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipse' must have a parameter 'Point' specified of format e.g.: '1.2,3.4'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_mat( tPoint, tPointMat );
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
                    moris::string_to_mat( tSemiDiameters, tSemiDiametersMat );
                    MORIS_ERROR( tSemiDiametersMat.numel() == 2 || tSemiDiametersMat.n_cols() == 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'SemiDiameters' vector for a 'ellipse' needs to be two." );

                    // get the exponent, default to 2.0 if not specified
                    std::string tExponent = "2.0";
                    mXmlReader->get_from_buffer( "Exponent", tExponent, std::string( "2.0" ) );
                    if ( tExponent == "" ) { tExponent = "2.0"; };

                    // set the parameters in the GEN parameter list
                    tGenParamList( 1 )( iGeom ).set( "type", "superellipse" );
                    std::string tParamString = tPoint + "," + tSemiDiameters + "," + tExponent + ",1.0,0.0,0.0";
                    tGenParamList( 1 )( iGeom ).set( "constant_parameters", tParamString );
                }

                // -------------------------------- //
                // SPHERE

                else if ( tPreDefGeom == "sphere" )
                {
                    // check dimensionality
                    MORIS_ERROR( tNumDims != 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "A 'sphere' can only be defined in 3D, but mesh is 2D." );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'circle' must have a parameter 'Point' specified of format e.g.: '1.2,3.4,5.6'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_mat( tPoint, tPointMat );
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
                    tGenParamList( 1 )( iGeom ).set( "type", "sphere" );
                    tGenParamList( 1 )( iGeom ).set( "constant_parameters", tPoint + "," + tRadius );
                }

                // -------------------------------- //
                // ELLIPSE

                else if ( tPreDefGeom == "ellipsoid" )
                {
                    // check dimensionality
                    MORIS_ERROR( tNumDims != 2,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "A 'sphere' can only be defined in 3D, but mesh is 2D." );

                    // get the point
                    std::string tPoint = "";
                    mXmlReader->get_from_buffer( "Point", tPoint, std::string( "" ) );
                    MORIS_ERROR( tPoint != "",
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "The 'ellipsoid' must have a parameter 'Point' specified of format e.g.: '1.2,3.4,5.6'" );
                    Matrix< DDRMat > tPointMat;
                    moris::string_to_mat( tPoint, tPointMat );
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
                    moris::string_to_mat( tSemiDiameters, tSemiDiametersMat );
                    MORIS_ERROR( tSemiDiametersMat.numel() == 3 || tSemiDiametersMat.n_cols() == 3,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Number of entries in 'SemiDiameters' vector for a 'ellipsoid' needs to be three." );

                    // get the exponent, default to 2.0 if not specified
                    std::string tExponent = "2.0";
                    mXmlReader->get_from_buffer( "Exponent", tExponent, std::string( "2.0" ) );
                    if ( tExponent == "" ) { tExponent = "2.0"; };

                    // set the parameters in the GEN parameter list
                    tGenParamList( 1 )( iGeom ).set( "type", "superellipsoid" );
                    tGenParamList( 1 )( iGeom ).set( "constant_parameters", tPoint + "," + tSemiDiameters + "," + tExponent );
                }

                // -------------------------------- //
                // UNKNOWN GEOMETRY

                else
                {
                    MORIS_ERROR( false,
                            "Library_IO_Meshgen::load_parameters_from_xml() - "
                            "Unknown pre-defined geometry. Supported Options are 'plane', 'circle', 'sphere', 'ellipse', and 'ellipsoid'." );
                }
            }    // end if: pre-defined geometry

            else if ( tGeomType == "image_file" )
            {
                // initialize with the image sdf default parameter list
                tGenParamList( 1 )( iGeom ) = prm::create_image_sdf_field_parameter_list();

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
                Matrix< DDRMat > tImageDimVec;
                moris::string_to_mat( tImageDimStr, tImageDimVec );
                MORIS_ERROR( tImageDimVec.numel() == tNumDims || tImageDimVec.n_cols() == tNumDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ImageDimensions' vector does not match number of spatial dimensions" );
                tGenParamList( 1 )( iGeom ).set( "image_dimensions", tImageDimStr );

                // get the image offset
                std::string tImageOffsetStr = "";
                mXmlReader->get_from_buffer( "ImageOrigin", tImageOffsetStr, std::string( "" ) );
                MORIS_ERROR( tFileName != "",
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "No image origin provided for geometry #%i. "
                        "Please provide image origin/offset using the tag 'ImageOrigin' in the the format e.g. '1.2,3.3' ",
                        iGeom );
                Matrix< DDRMat > tImageOffsetVec;
                moris::string_to_mat( tImageOffsetStr, tImageOffsetVec );
                MORIS_ERROR( tImageOffsetVec.numel() == tNumDims || tImageOffsetVec.n_cols() == tNumDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ImageOrigin' vector does not match number of spatial dimensions" );
                tGenParamList( 1 )( iGeom ).set( "image_offset", tImageOffsetStr );

            }    // end if: geometry from image file

            else if ( tGeomType == "object_file" )
            {
                // initialize with the sdf field default parameter list
                tGenParamList( 1 )( iGeom ) = prm::create_sdf_field_parameter_list();

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
                Matrix< DDRMat > tObjectOffsetVec;
                moris::string_to_mat( tObjectOffsetStr, tObjectOffsetVec );
                MORIS_ERROR( tObjectOffsetVec.numel() == tNumDims || tObjectOffsetVec.n_cols() == tNumDims,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Number of entries in 'ObjectOrigin' vector for geometry %i does not match number of spatial dimensions.",
                        iGeom );
                tGenParamList( 1 )( iGeom ).set( "sdf_object_offset", tObjectOffsetStr );

                // get an offset in the sign distance value if input is provided
                double tSdfShift = 0.0;
                mXmlReader->get_from_buffer( "SdfShift", tSdfShift, double( 0.0 ) );
                tGenParamList( 1 )( iGeom ).set( "sdf_shift", tSdfShift );
                
            }    // end if: geometry from image file

            else
            {
                MORIS_ERROR( false,
                        "Library_IO_Meshgen::load_parameters_from_xml() - "
                        "Geometry type '%s' unknown. Currently supported options are: 'pre_defined', 'image_file' and 'object_file'.",
                        tGeomType );
            }

        }    // end for: loop over geometries in input file

        // ------------------------------
        // Phase map

        // get the string from the input file
        std::string tPhaseMapString = "";
        mXmlReader->get( tGenPath + ".PhaseMap", tPhaseMapString, std::string( "" ) );
        bool tPhaseMapSpecified = ( tPhaseMapString != "" );

        // if a phase map has been specified, convert it to a matrix
        Matrix< DDUMat > tPhaseMap;
        moris::string_to_mat( tPhaseMapString, tPhaseMap );

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
    }

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_all_standard_parameters()
    {
        // go over all modules and check that
        for ( uint iModule = 0; iModule < (uint)( Parameter_List_Type::END_ENUM ); iModule++ )
        {
            // get the current module
            Parameter_List_Type tParamListType = (Parameter_List_Type)( iModule );

            // fill the parameter list entry with the standard parameters
            this->create_standard_parameter_list_for_module( tParamListType, mParameterLists( iModule ) );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_parameter_list_for_module(
            Parameter_List_Type  aParamListType,
            ModuleParameterList& aParameterList )
    {
        switch ( aParamListType )
        {
            case Parameter_List_Type::OPT:
                this->create_standard_OPT_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::HMR:
                this->create_standard_HMR_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::XTK:
                this->create_standard_XTK_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::GEN:
                this->create_standard_GEN_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::END_ENUM:
                MORIS_ERROR( false,
                        "Library_IO_Meshgen::create_standard_parameter_list_for_module() - "
                        "No standard library defined for module END_ENUM" );
                break;

            // create an empty parameter list for modules that are not needed
            default:
                aParameterList = ModuleParameterList();
                break;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_OPT_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_opt_problem_parameter_list();    // ParameterList();
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_XTK_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_xtk_parameter_list();    // ParameterList();

        // enrichment
        aParameterList( 0 )( 0 ).set( "enrich", true );
        aParameterList( 0 )( 0 ).set( "basis_rank", "bspline" );

        // output
        aParameterList( 0 )( 0 ).set( "output_file", "foreground_mesh.exo" );
        aParameterList( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterList( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterList( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterList( 0 )( 0 ).set( "global_T_matrix_output_file", "" );
        aParameterList( 0 )( 0 ).set( "elemental_T_matrix_output_file", "Elemental_Extraction_Operators" );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_HMR_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_hmr_parameter_list();    // ParameterList();

        // reduce buffer size as much as possible
        aParameterList( 0 )( 0 ).set( "refinement_buffer", 0 );
        aParameterList( 0 )( 0 ).set( "staircase_buffer", 0 );

        // Lagrange mesh is always 0
        aParameterList( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_GEN_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 3 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_gen_parameter_list();    // ParameterList();
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
