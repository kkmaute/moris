/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.cpp
 *
 */

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "fn_dot.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Parameters::Level_Set_Parameters( const Parameter_List& aParameterList )
            : Field_Parameters( aParameterList )
            , Design_Parameters( aParameterList )
            , mIsocontourThreshold( aParameterList.get< real >( "isocontour_threshold" ) )
            , mIsocontourTolerance( aParameterList.get< real >( "isocontour_tolerance" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Geometry::Level_Set_Geometry(
            std::shared_ptr< Field > aField,
            Level_Set_Parameters     aParameters,
            Node_Manager&            aNodeManager )
            : Geometry( aParameters, aParameters.mIntersectionTolerance )
            , Design_Field( aField, aParameters, aNodeManager )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = &aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_isocontour_threshold() const
    {
        return mParameters.mIsocontourThreshold;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Level_Set_Geometry::depends_on_advs() const
    {
        return mField->has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // If it's a background node, can get geometric region from field value
        if ( mNodeManager->is_background_node( aNodeIndex ) )
        {
            return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
        }
        else
        {
            // Get derived node
            const Node_Manager& tNodeManager = *mNodeManager;
            const Derived_Node& tDerivedNode = tNodeManager.get_derived_node( aNodeIndex );

            // If derived node knows it is on this interface, can return interface
            if ( tDerivedNode.is_on_interface( *this ) )
            {
                return Geometric_Region::INTERFACE;
            }
            else
            {
                // Get locators
                const Vector< Basis_Node >& tLocators = tDerivedNode.get_locator_nodes();

                // If we only have 2 locators, can use special logic if at least one node is on the interface
                if ( tLocators.size() == 2 )
                {
                    if ( this->get_geometric_region( tLocators( 0 ).get_index(), tLocators( 0 ).get_global_coordinates() ) == Geometric_Region::INTERFACE )
                    {
                        // Take geometric region of second node
                        return this->get_geometric_region( tLocators( 1 ).get_index(), tLocators( 1 ).get_global_coordinates() );
                    }
                    else if ( this->get_geometric_region( tLocators( 1 ).get_index(), tLocators( 1 ).get_global_coordinates() ) == Geometric_Region::INTERFACE )
                    {
                        // Take geometric region of first node
                        return this->get_geometric_region( tLocators( 0 ).get_index(), tLocators( 0 ).get_global_coordinates() );
                    }
                    else
                    {
                        // Resort to field value
                        return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
                    }
                }
                else
                {
                    // Must use field value for more than 2 locators
                    return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Level_Set_Geometry::create_intersection_node(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder )
    {
        if ( this->use_multilinear_interpolation() )
        {
            // Create multilinear intersection node
            return new Intersection_Node_Bilinear(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    *this );
        }
        else
        {
            // Create linear intersection node
            return new Intersection_Node_Linear(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    *this );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Level_Set_Geometry::compute_intersection_local_coordinate(
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode )
    {
        if ( this->use_multilinear_interpolation() )
        {
            // get isocontour threshold from geometry
            real tIsocontourThreshold = this->get_isocontour_threshold();

            // spatial dimension
            uint tNumDims = aFirstParentNode.get_global_coordinates().length();

            // number of nodes to be used for interpolation
            uint tNumBases;

            // build interpolator
            mtk::Interpolation_Function_Factory tFactory;
            mtk::Interpolation_Function_Base*   tInterpolation;

            // create interpolation function based on spatial dimension of problem
            switch ( tNumDims )
            {
                case 2:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 4;
                    break;
                }
                case 3:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 8;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Intersection_Node_Bilinear::compute_intersection - Interpolation type not implemented." );
                }
            }

            // allocate matrix for level set values at background cell nodes
            Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

            // get level set values of corner nodes
            for ( uint iBackgroundNode = 0; iBackgroundNode < tNumBases; ++iBackgroundNode )
            {
                tPhiBCNodes( iBackgroundNode ) = this->get_field_value(
                        aBackgroundNodes( iBackgroundNode )->get_index(),
                        aBackgroundNodes( iBackgroundNode )->get_global_coordinates() );
            }

            // Scale element level set field such that norm equals 1.0
            real tPhiScaling = 1.0 / norm( tPhiBCNodes );
            tPhiBCNodes      = tPhiScaling * tPhiBCNodes;
            tIsocontourThreshold *= tPhiScaling;

            // Get scaled parent level set values
            real tFirstParentPhi  = tPhiScaling * this->get_field_value( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() );
            real tSecondParentPhi = tPhiScaling * this->get_field_value( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() );

            // check that line is intersected
            if ( ( tFirstParentPhi - tIsocontourThreshold ) * ( tSecondParentPhi - tIsocontourThreshold ) > 0 )
            {
                delete tInterpolation;

                return MORIS_REAL_MAX;
            }

            // check that difference between first and second parent is not smaller than MORIS_REAL_EPS
            if ( std::abs( tFirstParentPhi - tSecondParentPhi ) < MORIS_REAL_EPS )
            {
                delete tInterpolation;

                // return that intersection node is at center of edge
                return 0.0;
            }

            // set Newton parameters
            const uint tNewMaxIter  = 20;    // maximum number of iterations in Newton
            const uint tCurvMaxIter = 10;    // maximum number of iterations using curvature information

            const real tNewRelax  = 1.0;      // relaxation factor for solution update
            const real tNewRelEps = 1e-8;     // required relative residual drop
            const real tNewAbsEps = 1e-12;    // required absolute residual drop
            const real tCurvMin   = 1e-8;     // threshold for curvature magnitude above which curvature is used

            // allocate matrices used within Newton loop
            Matrix< DDRMat > tCellCoordinate( tNumDims, 1 );
            Matrix< DDRMat > tBasis;
            Matrix< DDRMat > tDBasisDxi;
            Matrix< DDRMat > tDBasisDxi2;
            Matrix< DDRMat > tD2PhiDxi2;

            // vector from first to second parent in parametric coordinates
            Matrix< DDRMat > tSecondToFirstParent = aSecondParentNode.get_parametric_coordinates() - aFirstParentNode.get_parametric_coordinates();

            // initialized reference residual
            real tReferenceResidual = 0.0;
            real tResidual          = 0.0;

            // compute initial guess: location of intersection point along edge in edge CS
            real tEdgeCoordinate = ( 2.0 * tIsocontourThreshold - tFirstParentPhi - tSecondParentPhi )
                                 / ( tSecondParentPhi - tFirstParentPhi );
            Matrix< DDRMat > tInitialGuess = { { std::min( 1.0, std::max( tEdgeCoordinate, -1.0 ) ), -1.0, 1.0 } };

            // loop over initial guess trials
            for ( uint iGuess = 0; iGuess < tInitialGuess.numel(); iGuess++ )
            {
                // set initial guess
                tEdgeCoordinate = tInitialGuess( iGuess );

                // perform iterations
                for ( uint iNew = 0; iNew < tNewMaxIter; ++iNew )
                {
                    // compute local coordinate in background cell CS
                    tCellCoordinate = 0.5 * ( 1.0 - tEdgeCoordinate ) * aFirstParentNode.get_parametric_coordinates()
                                    + 0.5 * ( 1.0 + tEdgeCoordinate ) * aSecondParentNode.get_parametric_coordinates();

                    // compute basis function
                    tInterpolation->eval_N( tCellCoordinate, tBasis );

                    // compute residual
                    tResidual = dot( tBasis, tPhiBCNodes ) - tIsocontourThreshold;

                    // check convergence against absolute residual
                    if ( std::abs( tResidual ) < tNewAbsEps )
                    {
                        delete tInterpolation;

                        return tEdgeCoordinate;
                    }

                    // store reference residual
                    if ( iNew == 0 )
                    {
                        tReferenceResidual = std::abs( tResidual );
                    }

                    // check convergence against relative residual
                    if ( std::abs( tResidual ) < tNewRelEps * tReferenceResidual )
                    {
                        delete tInterpolation;

                        return tEdgeCoordinate;
                    }

                    // compute Jacobian
                    tInterpolation->eval_dNdXi( tCellCoordinate, tDBasisDxi );

                    // compute first order gradient of residual with respect to edge coordinate
                    real tGradRes = 0.5 * dot( tSecondToFirstParent, tDBasisDxi * tPhiBCNodes );

                    // initialize solution increment
                    real tSolIncrement;

                    // initialize curvature variable
                    real tSqrt2   = -1.0;
                    real tCurvRes = 0.0;

                    // compute second order derivatives
                    if ( iNew < tCurvMaxIter )
                    {
                        // compute Hessian
                        tInterpolation->eval_d2NdXi2( tCellCoordinate, tDBasisDxi2 );

                        tD2PhiDxi2 = tDBasisDxi2 * tPhiBCNodes;

                        // compute second order derivative of residual with respect to edge coordinate
                        if ( tNumDims == 2 )
                        {
                            tCurvRes =                                                                                   //
                                    0.125 * tD2PhiDxi2( 0 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 0 ) +    //
                                    0.125 * tD2PhiDxi2( 1 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 1 ) +    //
                                    0.250 * tD2PhiDxi2( 2 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 1 );
                        }
                        else
                        {
                            tCurvRes =                                                                                   //
                                    0.125 * tD2PhiDxi2( 0 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 0 ) +    //
                                    0.125 * tD2PhiDxi2( 1 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 1 ) +    //
                                    0.125 * tD2PhiDxi2( 2 ) * tSecondToFirstParent( 2 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 3 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 4 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 5 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 1 );
                        }

                        // compute square of sqrt term in root finding formula for quadratic equations
                        tSqrt2 = tGradRes * tGradRes - 4.0 * tCurvRes * tResidual;
                    }

                    // update solution depending on curvature
                    if ( iNew < tCurvMaxIter && std::abs( tCurvRes ) > tCurvMin && tSqrt2 >= 0.0 )
                    {
                        // add roots of quadratic function to current edge coordinate
                        real tSqrt = std::sqrt( tSqrt2 );

                        real tRootMin = ( -1.0 * tGradRes - tSqrt ) / ( 2.0 * tCurvRes );
                        real tRootMax = ( -1.0 * tGradRes + tSqrt ) / ( 2.0 * tCurvRes );

                        // check with root is in [-1,1] interval
                        bool tRootMinIsValid = std::abs( tEdgeCoordinate + tRootMin ) <= 1.0;
                        bool tRootMaxIsValid = std::abs( tEdgeCoordinate + tRootMax ) <= 1.0;

                        // compute increment
                        if ( tRootMinIsValid && tRootMaxIsValid )
                        {
                            tSolIncrement = std::abs( tRootMin ) < std::abs( tRootMax ) ? tRootMin : tRootMax;
                        }
                        else if ( tRootMinIsValid )
                        {
                            tSolIncrement = tRootMin;
                        }
                        else if ( tRootMaxIsValid )
                        {
                            tSolIncrement = tRootMax;
                        }
                        else
                        {
                            tSolIncrement = -1.0 * tResidual / ( tGradRes + MORIS_REAL_EPS );
                        }
                    }
                    else
                    {
                        // compute increment
                        tSolIncrement = -1.0 * tResidual / ( tGradRes + MORIS_REAL_EPS );
                    }

                    // update solution
                    tEdgeCoordinate += tNewRelax * tSolIncrement;

                    // trim solution
                    tEdgeCoordinate = std::min( 1.0, std::max( tEdgeCoordinate, -1.0 ) );
                }
            }

            // print debug information
            // for ( uint in = 0; in < tNumBases; ++in )
            // {
            //     std::string tStrg = "Anchestor_Node_" + std::to_string( aAncestorNodeIndices( in ) );
            //     print( aAncestorNodeCoordinates( in ), tStrg );
            // }

            // print( tPhiBCNodes, "tPhiBCNodes" );

            // fprintf( stderr, "tFirstParentPhi =%e   tSecondParentPhi = %e\n", tFirstParentPhi, tSecondParentPhi );

            // print( aFirstParentNodeLocalCoordinates, "aFirstParentNodeLocalCoordinates" );
            // print( aSecondParentNodeLocalCoordinates, "aSecondParentNodeLocalCoordinates" );

            MORIS_ERROR( false,
                    "Intersection_Node_Bilinear::compute_intersection - Newton did not converge: %s %e %s %e %s %e",
                    "Reference residual",
                    tReferenceResidual,
                    "Current residual",
                    std::abs( tResidual ),
                    "Edge coordinate",
                    tEdgeCoordinate );

            delete tInterpolation;

            return tEdgeCoordinate;
        }
        else
        {
            // Interface geometry values
            Matrix< DDRMat > tInterfaceGeometryValues = { { this->get_field_value( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() ) },
                { this->get_field_value( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() ) } };

            // Get isocontour threshold
            real tIsocontourThreshold = this->get_isocontour_threshold();

            // Interpolate
            Matrix< DDRMat > tLocalCoordinates = Interpolation::linear_interpolation_value( tInterfaceGeometryValues, tIsocontourThreshold );

            return tLocalCoordinates( 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::get_dfield_dcoordinates(
            const Basis_Node& aParentNode,
            Matrix< DDRMat >& aSensitivities ) const
    {
        if ( this->use_multilinear_interpolation() )
        {
            // TODO
            MORIS_ERROR( false, "Multilinear interpolation was detected for an intersection on another intersection. This needs to be implemented." );
        }
        else
        {
            // Get parents
            const Vector< Basis_Node >& tParentNodes = aParentNode.get_locator_nodes();

            // Geometry values
            real tDeltaPhi = this->get_field_value( tParentNodes( 1 ).get_index(), tParentNodes( 1 ).get_global_coordinates() )
                           - this->get_field_value( tParentNodes( 0 ).get_index(), tParentNodes( 0 ).get_global_coordinates() );

            // Compute parent vector
            Matrix< DDRMat > tParentVector = tParentNodes( 1 ).get_global_coordinates() - tParentNodes( 0 ).get_global_coordinates();

            // get number of spatial dimensions;
            uint tNumDim = tParentVector.length();

            // Compute square of length of parent vector
            real tParentLengthSquared = dot( tParentVector, tParentVector );

            // Sensitivities: dPhi/dx_i  = delta(Phi) / L_i where L_i = PaerentVectorLenth^2 / (ParentVector * e_i)
            for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumDim; tCoordinateIndex++ )
            {
                aSensitivities( tCoordinateIndex ) = tDeltaPhi * tParentVector( tCoordinateIndex ) / ( tParentLengthSquared + MORIS_REAL_EPS );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > > Level_Set_Geometry::get_mtk_fields()
    {
        // Get MTK field
        std::shared_ptr< mtk::Field > tMTKField = this->get_mtk_field();

        // Add to cell if it exists
        if ( tMTKField )
        {
            return { this->get_mtk_field() };
        }
        else
        {
            return {};
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::determine_geometric_region( real aLevelSetValue ) const
    {
        // Determine if value indicates that this point is on the interface
        if ( std::abs( aLevelSetValue - mParameters.mIsocontourThreshold ) <= mParameters.mIsocontourTolerance )
        {
            return Geometric_Region::INTERFACE;
        }

        // Otherwise, give the region relative to the isocontour threshold
        else if ( aLevelSetValue - mParameters.mIsocontourThreshold < 0 )
        {
            return Geometric_Region::NEGATIVE;
        }
        else
        {
            return Geometric_Region::POSITIVE;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Level_Set_Geometry::get_name()
    {
        return Design_Field::get_name();
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::string >
    Level_Set_Geometry::get_field_names()
    {
        Vector< std::string > tFieldNames( 1 );

        // Assign a default if this geometry does not have a name
        this->get_name() == "" ? tFieldNames( 0 ) = "Level Set Geometry: " + std::to_string( mOffsetID ) : tFieldNames( 0 ) = this->get_name();

        return tFieldNames;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        return Design_Field::import_advs( aOwnedADVs );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh )
    {
        Design_Field::reset_nodal_data( aInterpolationMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::discretize(
            mtk::Mesh_Pair    aMeshPair,
            sol::Dist_Vector* aOwnedADVs )
    {
        std::cout << aMeshPair.get_interpolation_mesh()->get_num_nodes(); // brendan

        if ( mSharedADVIDs.size() == 0 )
        {
            Design_Field::discretize( aMeshPair, aOwnedADVs, { {} }, mOffsetID );
        }
        else
        {
            MORIS_ASSERT( mSharedADVIDs.size() == 1,
                    "discretize() - Level Set geometries should have max one set of shared ADV IDs. Size = %ld",
                    mSharedADVIDs.size() );
            Design_Field::discretize( aMeshPair, aOwnedADVs, mSharedADVIDs( 0 ), mOffsetID );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Level_Set_Geometry::discretize(
            std::shared_ptr< mtk::Field > aMTKField,
            mtk::Mesh_Pair                aMeshPair,
            sol::Dist_Vector*             aOwnedADVs )
    {
        std::cout << aMeshPair.get_interpolation_mesh()->get_num_nodes(); // brendan

        if ( aMTKField->get_label() == this->get_name() )
        {
            if ( mSharedADVIDs.size() == 0 )
            {
                Design_Field::discretize( aMTKField, aMeshPair, aOwnedADVs, { {} }, mOffsetID );
            }
            else
            {
                MORIS_ASSERT( mSharedADVIDs.size() == 1,
                        "discretize() - Level Set geometries should have max one set of shared ADV IDs. Size = %ld",
                        mSharedADVIDs.size() );
                Design_Field::discretize( aMTKField, aMeshPair, aOwnedADVs, mSharedADVIDs( 0 ), mOffsetID );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::get_design_info(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Vector< real >&         aOutputDesignInfo )
    {
        aOutputDesignInfo.resize( 1 );
        aOutputDesignInfo( 0 ) = Design_Field::get_field_value( aNodeIndex, aCoordinates );
    }

    bool Level_Set_Geometry::intended_discretization()
    {
        return ( mParameters.mDiscretizationIndex >= 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Level_Set_Geometry::get_discretization_mesh_index()
    {
        MORIS_ASSERT( mParameters.mDiscretizationIndex >= 0,
                "A discretization is not intended for this field. Check this with intended_discretization() first." );

        return mParameters.mDiscretizationIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_discretization_lower_bound()
    {
        return mParameters.mDiscretizationLowerBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_discretization_upper_bound()
    {
        return mParameters.mDiscretizationUpperBound;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::update_dependencies( Vector< std::shared_ptr< Design > > aAllUpdatedDesigns )
    {
        // Get fields from designs
        Vector< std::shared_ptr< Field > > tUpdatedFields( aAllUpdatedDesigns.size() );
        for ( uint iFieldIndex = 0; iFieldIndex < tUpdatedFields.size(); iFieldIndex++ )
        {
            tUpdatedFields.append( aAllUpdatedDesigns( iFieldIndex )->get_fields() );
        }

        // Update fields
        Design_Field::update_dependencies( tUpdatedFields );
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
