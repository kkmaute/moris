
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp" //FEM/MSI/src

#define protected public
#define private   public
#include "cl_FEM_Set.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp" //FEM/INT/src
#include "cl_FEM_Element.hpp" //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

        // This test checks that a fem::Element computes its volume and size properly.
        TEST_CASE( "FEM Element Volume", "[moris],[fem],[FEM_Element_Volume]" )
        {
            // data for the fem::Element
            fem::Set     tFEMSet;

            // create a geometry interpolator for the FEM Set
            // define an integration mesh
            //------------------------------------------------------------------------------
            // integration mesh geometry type
            mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::QUAD;

            // define a QUAD4 integration element, i.e. space param coordinates xiHat
            Matrix< DDRMat > tXiHatIG  = {{ -1.0, -1.0 }, {  0.0, -1.0 }, {  0.0,  1.0 }, { -1.0,  1.0 }};
            Matrix< DDRMat > tTauHatIG = {{-1.0}, {1.0}, {0.0}};

            // the QUAD4 integration element in space physical coordinates xHat
            Matrix< DDRMat > tXHatIG = {{ 0.0, 0.0 }, { 0.5, 0.0 }, { 0.5, 1.0 }, { 0.0, 1.0 }};
            Matrix< DDRMat > tTHatIG = {{ 0.0 }, { 1.0 }, { 0.5 }};

            // integration mesh interpolation rule
            Interpolation_Rule tGeoInterpIGRule( tGeoTypeIG,
                                                 Interpolation_Type::LAGRANGE,
                                                 mtk::Interpolation_Order::LINEAR,
                                                 Interpolation_Type::LAGRANGE,
                                                 mtk::Interpolation_Order::QUADRATIC );

            // create a space and time geometry interpolator fot the integration element
            Geometry_Interpolator* tGeoInterpIG = new Geometry_Interpolator( tGeoInterpIGRule );

            //set the coefficients xHat, tHat
            tGeoInterpIG->set_space_coeff( tXHatIG );
            tGeoInterpIG->set_time_coeff(  tTHatIG );

            //set the coefficients xiHat, tauHat
            tGeoInterpIG->set_space_param_coeff( tXiHatIG );
            tGeoInterpIG->set_time_param_coeff(  tTauHatIG );

            tFEMSet.mMasterIGGeometryInterpolator = tGeoInterpIG;

            // create a integration rule
            Integration_Rule tIntegrationRule( tGeoTypeIG,
                                               Integration_Type::GAUSS,
                                               Integration_Order::QUAD_3x3,
                                               Integration_Type::GAUSS,
                                               Integration_Order::BAR_3 );

            // create a side integrator
            Integrator tIntegrator( tIntegrationRule );

            //get number of integration points, integration points and weights
            uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );
            Matrix< DDRMat > tIntegWeights;
            tIntegrator.get_weights( tIntegWeights );

            tFEMSet.mIntegPoints  = tIntegPoints;
            tFEMSet.mIntegWeights = tIntegWeights;

            // create a fem::Element
            fem::Element_Bulk tFEMElement;

            // set the element set
            tFEMElement.mSet = &tFEMSet;

            // compute element volume and check
            real tVolume = tFEMElement.compute_volume();
            CHECK( tVolume - 0.5 < 1E-6 );

            // compute element size and check
            real tSize = tFEMElement.compute_size( 2 );
            CHECK( tSize - 0.797885 < 1E-6 );

        }/* TEST_CASE */

        // This test checks that a fem::Cluster computes its volume properly.
        TEST_CASE( "FEM Cluster Volume", "[moris],[fem],[FEM_Cluster_Volume]" )
        {
            // data for the fem::Element
            fem::Set     tFEMSet;

            // create a geometry interpolator for the FEM Set
            // integration mesh geometry type
            mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::QUAD;

            // define a QUAD4 integration element, i.e. space param coordinates xiHat
            Matrix< DDRMat > tXiHatIG  = {{ -1.0, -1.0 }, {  0.0, -1.0 }, {  0.0,  1.0 }, { -1.0,  1.0 }};
            Matrix< DDRMat > tTauHatIG = {{-1.0}, {1.0}, {0.0}};

            // the QUAD4 integration element in space physical coordinates xHat
            Matrix< DDRMat > tXHatIG = {{ 0.0, 0.0 }, { 0.5, 0.0 }, { 0.5, 1.0 }, { 0.0, 1.0 }};
            Matrix< DDRMat > tTHatIG = {{ 0.0 }, { 1.0 }, { 0.5 }};

            // integration mesh interpolation rule
            Interpolation_Rule tGeoInterpIGRule( tGeoTypeIG,
                                                 Interpolation_Type::LAGRANGE,
                                                 mtk::Interpolation_Order::LINEAR,
                                                 Interpolation_Type::LAGRANGE,
                                                 mtk::Interpolation_Order::QUADRATIC );

            // create a space and time geometry interpolator fot the integration element
            Geometry_Interpolator* tGeoInterpIG = new Geometry_Interpolator( tGeoInterpIGRule );

            //set the coefficients xHat, tHat
            tGeoInterpIG->set_space_coeff( tXHatIG );
            tGeoInterpIG->set_time_coeff(  tTHatIG );

            //set the coefficients xiHat, tauHat
            tGeoInterpIG->set_space_param_coeff( tXiHatIG );
            tGeoInterpIG->set_time_param_coeff(  tTauHatIG );

            tFEMSet.mMasterIGGeometryInterpolator = tGeoInterpIG;

            // create a integration rule
            Integration_Rule tIntegrationRule( tGeoTypeIG,
                                               Integration_Type::GAUSS,
                                               Integration_Order::QUAD_3x3,
                                               Integration_Type::GAUSS,
                                               Integration_Order::BAR_3 );

            // create a side integrator
            Integrator tIntegrator( tIntegrationRule );

            //get number of integration points, integration points and weights
            uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );
            Matrix< DDRMat > tIntegWeights;
            tIntegrator.get_weights( tIntegWeights );

            // set the integration info in the set
            tFEMSet.mIntegPoints  = tIntegPoints;
            tFEMSet.mIntegWeights = tIntegWeights;

            // create a fem::Element
            Cell< fem::Element * > tFEMElements( 5, nullptr );

            // create the cluster elements
            for ( uint iElem = 0; iElem < 5; iElem++ )
            {
                // create a bulk element
                tFEMElements( iElem ) = new Element_Bulk();

                // set a fem set for the element
                tFEMElements( iElem )->mSet = &tFEMSet;
            }

            // create a cluster
            fem::Cluster tFEMCluster;

            // set the cluster elements
            tFEMCluster.mElements = tFEMElements;

            // compute cluster volume and check
            real tClusterVolume = tFEMCluster.compute_volume();
            CHECK( tClusterVolume - 2.5 < 1E-6 );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
