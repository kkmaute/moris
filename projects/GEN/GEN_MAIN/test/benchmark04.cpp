#include <string>
#include <iostream>
#include <cmath>
#include <utility>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_PRM_FEM_Parameters.hpp"
#include "cl_PRM_MSI_Parameters.hpp"
#include "cl_PRM_SOL_Parameters.hpp"
#include "cl_PRM_VIS_Parameters.hpp"
#include "cl_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "AztecOO.h"


#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /*
     * @brief apply traction to top and bottom, fix portion of right side away from crack
     * 
     *        crack - angled center crack with two tips (NOTE: known problems with decomposition when crack is at 45deg)
     */
    
void Func1( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}
//------------------------------------------------------------------------------
void Func2D( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
             moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
             moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = { { aParameters( 0 )( 0 ),                      0.0},
                    {                   0.0,    aParameters( 0 )( 1 )} };

}
//------------------------------------------------------------------------------
void Func2D_end( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                 moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                 moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix.set_size(2,2, 0.0);
    Matrix< DDRMat > tCoords = aFIManager->get_IP_geometry_interpolator()->valx();
    
    if( tCoords(1) > 4.5 && tCoords(1) < 5.5 )
    {
        aPropMatrix = { { aParameters( 0 )( 0 ),                      0.0},
                        {                   0.0,    aParameters( 0 )( 1 )} };
    }
}
//------------------------------------------------------------------------------
bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
{
    return true;
}
//------------------------------------------------------------------------------
void Lvl_set_0(       moris::real                & aReturnValue,
                const moris::Matrix< DDRMat >    & aPoint,
                const moris::Cell< moris::real > & aConst )
{
    real pi = 3.14159265358979323846;
    
    real betta = 30*(pi/180);   // angle in radians
    
    Matrix< DDRMat > tCenter_0(1, 2, 0.0);
//     tCenter_0(0) = 0.0;    tCenter_0(1) = 5.001;
    tCenter_0(0) = 5.001;    tCenter_0(1) = 5.001;
 
    Matrix< DDRMat > tNormal_0(1,2, 0.0);
//     tNormal_0(0) = 0.0;   tNormal_0(1) = 1.0;
    tNormal_0(0) = -std::sin(betta);   tNormal_0(1) = std::cos(betta);

    // plane which represents the crack
    moris::real tVal_0 = tNormal_0(0)*(aPoint(0)-tCenter_0(0)) + tNormal_0(1)*(aPoint(1)-tCenter_0(1));
    
    aReturnValue = tVal_0;
}
//------------------------------------------------------------------------------
void Lvl_set_1(       moris::real                & aReturnValue,
                const moris::Matrix< DDRMat >    & aPoint,
                const moris::Cell< moris::real > & aConst )
{  
    real pi = 3.14159265358979323846;
    
    real betta = 30*(pi/180);   // angle in radians
    
    Matrix< DDRMat > tCenter_1(1,2, 0.0);
    tCenter_1(0) = 5.501;    tCenter_1(1) = 5.001;
 
    Matrix< DDRMat > tNormal_1(1,2, 0.0);
//     tNormal_1(0) = 1.0;   tNormal_1(1) = 0.0;
    tNormal_1(0) = std::cos(betta);   tNormal_1(1) = std::sin(betta);

    // plane which represents the right tip of the crack
    moris::real tVal_1 = tNormal_1(0)*( aPoint(0)-tCenter_1(0) ) + tNormal_1(1)*( aPoint(1)-tCenter_1(1) );
    
    aReturnValue = tVal_1;
}
//------------------------------------------------------------------------------
void Lvl_set_2(       moris::real                & aReturnValue,
                const moris::Matrix< DDRMat >    & aPoint,
                const moris::Cell< moris::real > & aConst )
{
    real pi = 3.14159265358979323846;
    
    real betta = 30*(pi/180);   // angle in radians
    
    Matrix< DDRMat > tCenter_2(1,2, 0.0);
    tCenter_2(0) = 4.501;    tCenter_2(1) = 5.001;
    
    Matrix< DDRMat > tNormal_2(1,2, 0.0);
//     tNormal_2(0) = -1.0;   tNormal_2(1) = 0.0;
    tNormal_2(0) = std::cos(betta);   tNormal_2(1) = std::sin(betta);
    
    // plane which represents the left tip of the crack
    moris::real tVal_2 = tNormal_2(0)*( aPoint(0)-tCenter_2(0) ) + tNormal_2(1)*( aPoint(1)-tCenter_2(1) );   
    
    aReturnValue = tVal_2;
}
//------------------------------------------------------------------------------
void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
{
    // create a cell of cell of parameter list for fem
    tParameterList.resize( 5 );

    //------------------------------------------------------------------------------
    // fill the property part of the parameter list
    moris::uint tNumProperties = 16;
    tParameterList( 0 ).resize( tNumProperties );

    // create parameter list for property 1
    tParameterList( 0 )( 0 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 0 ).set( "property_name",            std::string("PropEModPlate0") );
    tParameterList( 0 )( 0 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 0 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 0 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 2
    tParameterList( 0 )( 1 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 1 ).set( "property_name",            std::string("PropEModPlate1") );
    tParameterList( 0 )( 1 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 1 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 1 ).set( "value_function",           std::string("Func1") );

    // create parameter list for property 3
    tParameterList( 0 )( 2 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 2 ).set( "property_name",            std::string("PropEModPlate3") );
    tParameterList( 0 )( 2 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 2 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 2 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 4
    tParameterList( 0 )( 3 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 3 ).set( "property_name",            std::string("PropEModPlate4") );
    tParameterList( 0 )( 3 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 3 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 3 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 5
    tParameterList( 0 )( 4 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 4 ).set( "property_name",            std::string("PropEModPlate5") );
    tParameterList( 0 )( 4 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 4 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 4 ).set( "value_function",           std::string("Func1") );    
    
    // create parameter list for property 6
    tParameterList( 0 )( 5 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 5 ).set( "property_name",            std::string("PropEModPlate7") );
    tParameterList( 0 )( 5 ).set( "function_parameters",      std::string("10.0") );
//     tParameterList( 0 )( 5 ).set( "function_parameters",      std::string("10000000.0") );
    tParameterList( 0 )( 5 ).set( "value_function",           std::string("Func1") );        
       
    //-----------------------------------------------
    
    // create parameter list for property 7
    tParameterList( 0 )( 6 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 6 ).set( "property_name",            std::string("PropNuPlate0") );
    tParameterList( 0 )( 6 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 6 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 6 ).set( "value_function",           std::string("Func1") );

    // create parameter list for property 8
    tParameterList( 0 )( 7 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 7 ).set( "property_name",            std::string("PropNuPlate1") );
    tParameterList( 0 )( 7 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 7 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 7 ).set( "value_function",           std::string("Func1") );

    // create parameter list for property 9
    tParameterList( 0 )( 8 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 8 ).set( "property_name",            std::string("PropNuPlate3") );
    tParameterList( 0 )( 8 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 8 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 8 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 10
    tParameterList( 0 )( 9 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 9 ).set( "property_name",            std::string("PropNuPlate4") );
    tParameterList( 0 )( 9 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 9 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 9 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 11
    tParameterList( 0 )( 10 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 10 ).set( "property_name",            std::string("PropNuPlate5") );
    tParameterList( 0 )( 10 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 10 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 10 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 12
    tParameterList( 0 )( 11 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 11 ).set( "property_name",            std::string("PropNuPlate7") );
    tParameterList( 0 )( 11 ).set( "function_parameters",      std::string("0.3") );
//     tParameterList( 0 )( 11 ).set( "function_parameters",      std::string("0.334") );
    tParameterList( 0 )( 11 ).set( "value_function",           std::string("Func1") );
    
    //-----------------------------------------------
    
    // create parameter list for property 13
    tParameterList( 0 )( 12 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 12 ).set( "property_name",            std::string("PropNeumannTop") );
    tParameterList( 0 )( 12 ).set( "function_parameters",      std::string("0.0;1.0") );
//     tParameterList( 0 )( 12 ).set( "function_parameters",      std::string("0.0;230000.0") );
    tParameterList( 0 )( 12 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 14
    tParameterList( 0 )( 13 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 13 ).set( "property_name",            std::string("PropNeumannBot") );
    tParameterList( 0 )( 13 ).set( "function_parameters",      std::string("0.0;-1.0") );
//     tParameterList( 0 )( 13 ).set( "function_parameters",      std::string("0.0;-230000.0") );
    tParameterList( 0 )( 13 ).set( "value_function",           std::string("Func1") );
    
    //-----------------------------------------------
    
    // create parameter list for property 15
    tParameterList( 0 )( 14 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 14 ).set( "property_name",            std::string("PropDirichletFixed") );
    tParameterList( 0 )( 14 ).set( "function_parameters",      std::string("0.0;0.0") );
    tParameterList( 0 )( 14 ).set( "value_function",           std::string("Func1") );
    
    // create parameter list for property 16
    tParameterList( 0 )( 15 ) = prm::create_property_parameter_list();
    tParameterList( 0 )( 15 ).set( "property_name",            std::string("PropDirichletFixedSelect") );
    tParameterList( 0 )( 15 ).set( "function_parameters",      std::string("1.0;1.0") );
//     tParameterList( 0 )( 15 ).set( "value_function",           std::string("Func2D") );
    tParameterList( 0 )( 15 ).set( "value_function",           std::string("Func2D_end") );
    
    //------------------------------------------------------------------------------
    // fill the constitutive model part of the parameter list
    moris::uint tNumCMs = 6;
    tParameterList( 1 ).resize( tNumCMs );

    // create parameter list for constitutive model 1
    tParameterList( 1 )( 0 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 0 ).set( "constitutive_name", std::string("CMPlate0") );
    tParameterList( 1 )( 0 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 0 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 0 ).set( "properties",        std::string("PropEModPlate0,YoungsModulus;PropNuPlate0,PoissonRatio") );
    tParameterList( 1 )( 0 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    // create parameter list for constitutive model 2
    tParameterList( 1 )( 1 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 1 ).set( "constitutive_name", std::string("CMPlate1") );
    tParameterList( 1 )( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 1 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 1 ).set( "properties",        std::string("PropEModPlate1,YoungsModulus;PropNuPlate1,PoissonRatio") );
    tParameterList( 1 )( 1 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    // create parameter list for constitutive model 3
    tParameterList( 1 )( 2 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 2 ).set( "constitutive_name", std::string("CMPlate3") );
    tParameterList( 1 )( 2 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 2 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 2 ).set( "properties",        std::string("PropEModPlate3,YoungsModulus;PropNuPlate3,PoissonRatio") );
    tParameterList( 1 )( 2 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    // create parameter list for constitutive model 4
    tParameterList( 1 )( 3 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 3 ).set( "constitutive_name", std::string("CMPlate4") );
    tParameterList( 1 )( 3 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 3 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 3 ).set( "properties",        std::string("PropEModPlate4,YoungsModulus;PropNuPlate4,PoissonRatio") );
    tParameterList( 1 )( 3 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    // create parameter list for constitutive model 5
    tParameterList( 1 )( 4 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 4 ).set( "constitutive_name", std::string("CMPlate5") );
    tParameterList( 1 )( 4 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 4 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 4 ).set( "properties",        std::string("PropEModPlate5,YoungsModulus;PropNuPlate5,PoissonRatio") );
    tParameterList( 1 )( 4 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    // create parameter list for constitutive model 6
    tParameterList( 1 )( 5 ) = prm::create_constitutive_model_parameter_list();
    tParameterList( 1 )( 5 ).set( "constitutive_name", std::string("CMPlate7") );
    tParameterList( 1 )( 5 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
    tParameterList( 1 )( 5 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
    tParameterList( 1 )( 5 ).set( "properties",        std::string("PropEModPlate7,YoungsModulus;PropNuPlate7,PoissonRatio") );
    tParameterList( 1 )( 5 ).set( "model_type",        static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
    
    //------------------------------------------------------------------------------
    // fill the stabilization parameter part of the parameter list
    moris::uint tNumSPs = 19;
    tParameterList( 2 ).resize( tNumSPs );

    // create parameter list for stabilization parameter 1
    tParameterList( 2 )( 0 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 0 ).set( "stabilization_name",  std::string("SPDirichletNitscheFixed_p3") );
    tParameterList( 2 )( 0 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
    tParameterList( 2 )( 0 ).set( "function_parameters", std::string("1000.0") );
    tParameterList( 2 )( 0 ).set( "master_properties",   std::string("PropEModPlate3,Material") );
    
   //-----------------------------------------------         
    
    tParameterList( 2 )( 1 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 1 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_0_1") );
    tParameterList( 2 )( 1 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 1 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 1 ).set( "master_properties",   std::string("PropEModPlate0,Material") );
    tParameterList( 2 )( 1 ).set( "slave_properties",    std::string("PropEModPlate1,Material") );
    
    tParameterList( 2 )( 2 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 2 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_0_1") );
    tParameterList( 2 )( 2 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 2 ).set( "master_properties",   std::string("PropEModPlate0,Material") );
    tParameterList( 2 )( 2 ).set( "slave_properties",    std::string("PropEModPlate1,Material") );
    
    tParameterList( 2 )( 3 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 3 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_0_1") );
    tParameterList( 2 )( 3 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 3 ).set( "master_properties",   std::string("PropEModPlate0,Material") );
    tParameterList( 2 )( 3 ).set( "slave_properties",    std::string("PropEModPlate1,Material") );
    
    //-----------------------------------------------         
    
    tParameterList( 2 )( 4 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 4 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_1_3") );
    tParameterList( 2 )( 4 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 4 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 4 ).set( "master_properties",   std::string("PropEModPlate1,Material") );
    tParameterList( 2 )( 4 ).set( "slave_properties",    std::string("PropEModPlate3,Material") );
    
    tParameterList( 2 )( 5 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 5 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_1_3") );
    tParameterList( 2 )( 5 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 5 ).set( "master_properties",   std::string("PropEModPlate1,Material") );
    tParameterList( 2 )( 5 ).set( "slave_properties",    std::string("PropEModPlate3,Material") );
    
    tParameterList( 2 )( 6 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 6 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_1_3") );
    tParameterList( 2 )( 6 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 6 ).set( "master_properties",   std::string("PropEModPlate1,Material") );
    tParameterList( 2 )( 6 ).set( "slave_properties",    std::string("PropEModPlate3,Material") );
    
    //-----------------------------------------------         
    
    tParameterList( 2 )( 7 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 7 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_3_7") );
    tParameterList( 2 )( 7 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 7 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 7 ).set( "master_properties",   std::string("PropEModPlate3,Material") );
    tParameterList( 2 )( 7 ).set( "slave_properties",    std::string("PropEModPlate7,Material") );
    
    tParameterList( 2 )( 8 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 8 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_3_7") );
    tParameterList( 2 )( 8 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 8 ).set( "master_properties",   std::string("PropEModPlate3,Material") );
    tParameterList( 2 )( 8 ).set( "slave_properties",    std::string("PropEModPlate7,Material") );

    tParameterList( 2 )( 9 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 9 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_3_7") );
    tParameterList( 2 )( 9 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 9 ).set( "master_properties",   std::string("PropEModPlate3,Material") );
    tParameterList( 2 )( 9 ).set( "slave_properties",    std::string("PropEModPlate7,Material") );
    
    //-----------------------------------------------         
    
    tParameterList( 2 )( 10 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 10 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_7_5") );
    tParameterList( 2 )( 10 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 10 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 10 ).set( "master_properties",   std::string("PropEModPlate7,Material") );
    tParameterList( 2 )( 10 ).set( "slave_properties",    std::string("PropEModPlate5,Material") );
    
    tParameterList( 2 )( 11 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 11 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_7_5") );
    tParameterList( 2 )( 11 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 11 ).set( "master_properties",   std::string("PropEModPlate7,Material") );
    tParameterList( 2 )( 11 ).set( "slave_properties",    std::string("PropEModPlate5,Material") );

    tParameterList( 2 )( 12 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 12 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_7_5") );
    tParameterList( 2 )( 12 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 12 ).set( "master_properties",   std::string("PropEModPlate7,Material") );
    tParameterList( 2 )( 12 ).set( "slave_properties",    std::string("PropEModPlate5,Material") );
    
    //-----------------------------------------------         
    
    tParameterList( 2 )( 13 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 13 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_5_4") );
    tParameterList( 2 )( 13 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 13 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 13 ).set( "master_properties",   std::string("PropEModPlate5,Material") );
    tParameterList( 2 )( 13 ).set( "slave_properties",    std::string("PropEModPlate4,Material") );
    
    tParameterList( 2 )( 14 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 14 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_5_4") );
    tParameterList( 2 )( 14 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 14 ).set( "master_properties",   std::string("PropEModPlate5,Material") );
    tParameterList( 2 )( 14 ).set( "slave_properties",    std::string("PropEModPlate4,Material") );

    tParameterList( 2 )( 15 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 15 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_5_4") );
    tParameterList( 2 )( 15 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 15 ).set( "master_properties",   std::string("PropEModPlate5,Material") );
    tParameterList( 2 )( 15 ).set( "slave_properties",    std::string("PropEModPlate4,Material") );
    
    //-----------------------------------------------         
    
    tParameterList( 2 )( 16 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 16 ).set( "stabilization_name",  std::string("SPNitscheInterfacePlateToPlate_4_0") );
    tParameterList( 2 )( 16 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
    tParameterList( 2 )( 16 ).set( "function_parameters", std::string("100.0") );
    tParameterList( 2 )( 16 ).set( "master_properties",   std::string("PropEModPlate4,Material") );
    tParameterList( 2 )( 16 ).set( "slave_properties",    std::string("PropEModPlate0,Material") );
    
    tParameterList( 2 )( 17 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 17 ).set( "stabilization_name",  std::string("SPMasterWeightPlateToPlate_4_0") );
    tParameterList( 2 )( 17 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 17 ).set( "master_properties",   std::string("PropEModPlate4,Material") );
    tParameterList( 2 )( 17 ).set( "slave_properties",    std::string("PropEModPlate0,Material") );

    tParameterList( 2 )( 18 ) = prm::create_stabilization_parameter_parameter_list();
    tParameterList( 2 )( 18 ).set( "stabilization_name",  std::string("SPSlaveWeightPlateToPlate_4_0") );
    tParameterList( 2 )( 18 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
    tParameterList( 2 )( 18 ).set( "master_properties",   std::string("PropEModPlate4,Material") );
    tParameterList( 2 )( 18 ).set( "slave_properties",    std::string("PropEModPlate0,Material") );
    
    //------------------------------------------------------------------------------
    // fill the IWG part of the parameter list
    moris::uint tNumIWGs = 19;
    tParameterList( 3 ).resize( tNumIWGs );

    // create parameter list for IWG 1
    tParameterList( 3 )( 0 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 0 ).set( "IWG_name",                   std::string("IWGPlateBulk0") );
    tParameterList( 3 )( 0 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 0 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 0 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 0 ).set( "master_constitutive_models", std::string("CMPlate0,ElastLinIso") );
    tParameterList( 3 )( 0 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p0,HMR_dummy_c_p0") );
    
    // create parameter list for IWG 2
    tParameterList( 3 )( 1 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 1 ).set( "IWG_name",                   std::string("IWGPlateBulk1") );
    tParameterList( 3 )( 1 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 1 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 1 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 1 ).set( "master_constitutive_models", std::string("CMPlate1,ElastLinIso") );
    tParameterList( 3 )( 1 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1,HMR_dummy_c_p1") );
    
    // create parameter list for IWG 3
    tParameterList( 3 )( 2 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 2 ).set( "IWG_name",                   std::string("IWGPlateBulk3") );
    tParameterList( 3 )( 2 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 2 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 2 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 2 ).set( "master_constitutive_models", std::string("CMPlate3,ElastLinIso") );
    tParameterList( 3 )( 2 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p3,HMR_dummy_c_p3") );
    
    // create parameter list for IWG 4
    tParameterList( 3 )( 3 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 3 ).set( "IWG_name",                   std::string("IWGPlateBulk4") );
    tParameterList( 3 )( 3 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 3 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 3 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 3 ).set( "master_constitutive_models", std::string("CMPlate4,ElastLinIso") );
    tParameterList( 3 )( 3 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p4,HMR_dummy_c_p4") );
    
    // create parameter list for IWG 5
    tParameterList( 3 )( 4 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 4 ).set( "IWG_name",                   std::string("IWGPlateBulk5") );
    tParameterList( 3 )( 4 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 4 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 4 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 4 ).set( "master_constitutive_models", std::string("CMPlate5,ElastLinIso") );
    tParameterList( 3 )( 4 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p5,HMR_dummy_c_p5") );
    
    // create parameter list for IWG 6
    tParameterList( 3 )( 5 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 5 ).set( "IWG_name",                   std::string("IWGPlateBulk7") );
    tParameterList( 3 )( 5 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
    tParameterList( 3 )( 5 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 5 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 5 ).set( "master_constitutive_models", std::string("CMPlate7,ElastLinIso") );
    tParameterList( 3 )( 5 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p7,HMR_dummy_c_p7") );
    
    //-----------------------------------------------
    
    // create parameter list for IWG 7
    tParameterList( 3 )( 6 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 6 ).set( "IWG_name",                   std::string("IWGNeumannTopP4") );
    tParameterList( 3 )( 6 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 6 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 6 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 6 ).set( "master_properties",          std::string("PropNeumannTop,Neumann") );
    tParameterList( 3 )( 6 ).set( "mesh_set_names",             std::string("SideSet_3_n_p4,SideSet_3_c_p4") );

    // create parameter list for IWG 8
    tParameterList( 3 )( 7 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 7 ).set( "IWG_name",                   std::string("IWGNeumannTopP5") );
    tParameterList( 3 )( 7 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 7 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 7 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 7 ).set( "master_properties",          std::string("PropNeumannTop,Neumann") );
    tParameterList( 3 )( 7 ).set( "mesh_set_names",             std::string("SideSet_3_n_p5,SideSet_3_c_p5") );
    
    // create parameter list for IWG 9
    tParameterList( 3 )( 8 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 8 ).set( "IWG_name",                   std::string("IWGNeumannTopP7") );
    tParameterList( 3 )( 8 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 8 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 8 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 8 ).set( "master_properties",          std::string("PropNeumannTop,Neumann") );
    tParameterList( 3 )( 8 ).set( "mesh_set_names",             std::string("SideSet_3_n_p7,SideSet_3_c_p7") );

    //-----------------------------------------------
    
    // create parameter list for IWG 10
    tParameterList( 3 )( 9 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 9 ).set( "IWG_name",                   std::string("IWGNeumannBotP0") );
    tParameterList( 3 )( 9 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 9 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 9 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 9 ).set( "master_properties",          std::string("PropNeumannBot,Neumann") );
    tParameterList( 3 )( 9 ).set( "mesh_set_names",             std::string("SideSet_1_n_p0,SideSet_1_c_p0") );
    
    // create parameter list for IWG 11
    tParameterList( 3 )( 10 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 10 ).set( "IWG_name",                   std::string("IWGNeumannBotP1") );
    tParameterList( 3 )( 10 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 10 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 10 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 10 ).set( "master_properties",          std::string("PropNeumannBot,Neumann") );
    tParameterList( 3 )( 10 ).set( "mesh_set_names",             std::string("SideSet_1_n_p1,SideSet_1_c_p1") );
    
    // create parameter list for IWG 12
    tParameterList( 3 )( 11 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 11 ).set( "IWG_name",                   std::string("IWGNeumannBotP3") );
    tParameterList( 3 )( 11 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
    tParameterList( 3 )( 11 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 11 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 11 ).set( "master_properties",          std::string("PropNeumannBot,Neumann") );
    tParameterList( 3 )( 11 ).set( "mesh_set_names",             std::string("SideSet_1_n_p3,SideSet_1_c_p3") );
    
    //-----------------------------------------------
    
    // create parameter list for IWG 13
    tParameterList( 3 )( 12 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 12 ).set( "IWG_name",                   std::string("IWGDirichletFixedSs2P3") );
    tParameterList( 3 )( 12 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
    tParameterList( 3 )( 12 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 12 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 12 ).set( "master_properties",          std::string("PropDirichletFixed,Dirichlet;PropDirichletFixedSelect,Select") );
    tParameterList( 3 )( 12 ).set( "master_constitutive_models", std::string("CMPlate1,ElastLinIso") );
    tParameterList( 3 )( 12 ).set( "stabilization_parameters",   std::string("SPDirichletNitscheFixed_p3,DirichletNitsche") );
    tParameterList( 3 )( 12 ).set( "mesh_set_names",             std::string("SideSet_2_n_p3") );
    
    //-----------------------------------------------
    
    // create parameter list for IWG 14
    tParameterList( 3 )( 13 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 13 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_0_1") );
    tParameterList( 3 )( 13 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 13 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 13 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 13 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_0_1,NitscheInterface;SPMasterWeightPlateToPlate_0_1,MasterWeightInterface;SPSlaveWeightPlateToPlate_0_1,SlaveWeightInterface") );
    tParameterList( 3 )( 13 ).set( "master_constitutive_models", std::string("CMPlate0,ElastLinIso") );
    tParameterList( 3 )( 13 ).set( "slave_constitutive_models",  std::string("CMPlate1,ElastLinIso") );
    tParameterList( 3 )( 13 ).set( "mesh_set_names",             std::string("dbl_iside_p0_0_p1_1") );
    
    // create parameter list for IWG 15
    tParameterList( 3 )( 14 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 14 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_1_3") );
    tParameterList( 3 )( 14 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 14 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 14 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 14 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_1_3,NitscheInterface;SPMasterWeightPlateToPlate_1_3,MasterWeightInterface;SPSlaveWeightPlateToPlate_1_3,SlaveWeightInterface") );
    tParameterList( 3 )( 14 ).set( "master_constitutive_models", std::string("CMPlate1,ElastLinIso") );
    tParameterList( 3 )( 14 ).set( "slave_constitutive_models",  std::string("CMPlate3,ElastLinIso") );
    tParameterList( 3 )( 14 ).set( "mesh_set_names",             std::string("dbl_iside_p0_1_p1_3") );
    
    // create parameter list for IWG 16
    tParameterList( 3 )( 15 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 15 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_3_7") );
    tParameterList( 3 )( 15 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 15 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 15 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 15 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_3_7,NitscheInterface;SPMasterWeightPlateToPlate_3_7,MasterWeightInterface;SPSlaveWeightPlateToPlate_3_7,SlaveWeightInterface") );
    tParameterList( 3 )( 15 ).set( "master_constitutive_models", std::string("CMPlate3,ElastLinIso") );
    tParameterList( 3 )( 15 ).set( "slave_constitutive_models",  std::string("CMPlate7,ElastLinIso") );
    tParameterList( 3 )( 15 ).set( "mesh_set_names",             std::string("dbl_iside_p0_3_p1_7") );
    
    // create parameter list for IWG 17
    tParameterList( 3 )( 16 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 16 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_7_5") );
    tParameterList( 3 )( 16 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 16 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 16 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 16 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_7_5,NitscheInterface;SPMasterWeightPlateToPlate_7_5,MasterWeightInterface;SPSlaveWeightPlateToPlate_7_5,SlaveWeightInterface") );
    tParameterList( 3 )( 16 ).set( "master_constitutive_models", std::string("CMPlate7,ElastLinIso") );
    tParameterList( 3 )( 16 ).set( "slave_constitutive_models",  std::string("CMPlate5,ElastLinIso") );
    tParameterList( 3 )( 16 ).set( "mesh_set_names",             std::string("dbl_iside_p0_5_p1_7") );
    
    // create parameter list for IWG 18
    tParameterList( 3 )( 17 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 17 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_5_4") );
    tParameterList( 3 )( 17 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 17 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 17 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 17 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_5_4,NitscheInterface;SPMasterWeightPlateToPlate_5_4,MasterWeightInterface;SPSlaveWeightPlateToPlate_5_4,SlaveWeightInterface") );
    tParameterList( 3 )( 17 ).set( "master_constitutive_models", std::string("CMPlate5,ElastLinIso") );
    tParameterList( 3 )( 17 ).set( "slave_constitutive_models",  std::string("CMPlate4,ElastLinIso") );
    tParameterList( 3 )( 17 ).set( "mesh_set_names",             std::string("dbl_iside_p0_4_p1_5") );
    
    // create parameter list for IWG 19
    tParameterList( 3 )( 18 ) = prm::create_IWG_parameter_list();
    tParameterList( 3 )( 18 ).set( "IWG_name",                   std::string("IWGInterfacePlateToPlate_4_0") );
    tParameterList( 3 )( 18 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE ) );
    tParameterList( 3 )( 18 ).set( "dof_residual",               std::string("UX,UY") );
    tParameterList( 3 )( 18 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 3 )( 18 ).set( "stabilization_parameters",   std::string("SPNitscheInterfacePlateToPlate_4_0,NitscheInterface;SPMasterWeightPlateToPlate_4_0,MasterWeightInterface;SPSlaveWeightPlateToPlate_4_0,SlaveWeightInterface") );
    tParameterList( 3 )( 18 ).set( "master_constitutive_models", std::string("CMPlate4,ElastLinIso") );
    tParameterList( 3 )( 18 ).set( "slave_constitutive_models",  std::string("CMPlate0,ElastLinIso") );
    tParameterList( 3 )( 18 ).set( "mesh_set_names",             std::string("dbl_iside_p0_0_p1_4") );
    
    //------------------------------------------------------------------------------
    // dbl_iside_p0_0_p1_2
    //------------------------------------------------------------------------------
    // fill the IQI part of the parameter list
    moris::uint tNumIQIs = 2;
    tParameterList( 4 ).resize( tNumIQIs );

    // create parameter list for IQI 1
    tParameterList( 4 )( 0 ) = prm::create_IQI_parameter_list();
    tParameterList( 4 )( 0 ).set( "IQI_name",                   std::string("IQIU_X") );
    tParameterList( 4 )( 0 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( 0 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UX ) );
    tParameterList( 4 )( 0 ).set( "vectorial_field_index",      0 );
    tParameterList( 4 )( 0 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 4 )( 0 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1,HMR_dummy_n_p3,HMR_dummy_c_p3,HMR_dummy_n_p4,HMR_dummy_c_p4,HMR_dummy_n_p5,HMR_dummy_c_p5,HMR_dummy_n_p7,HMR_dummy_c_p7") );
    
    // create parameter list for IQI 2
    tParameterList( 4 )( 1 ) = prm::create_IQI_parameter_list();
    tParameterList( 4 )( 1 ).set( "IQI_name",                   std::string("IQIU_Y") );
    tParameterList( 4 )( 1 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
    tParameterList( 4 )( 1 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UY ) );
    tParameterList( 4 )( 1 ).set( "vectorial_field_index",      1 );
    tParameterList( 4 )( 1 ).set( "master_dof_dependencies",    std::string("UX,UY") );
    tParameterList( 4 )( 1 ).set( "mesh_set_names",             std::string("HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1,HMR_dummy_n_p3,HMR_dummy_c_p3,HMR_dummy_n_p4,HMR_dummy_c_p4,HMR_dummy_n_p5,HMR_dummy_c_p5,HMR_dummy_n_p7,HMR_dummy_c_p7") );
    

//     tParameterList( 4 )( 0 ) = prm::create_IQI_parameter_list();
//     tParameterList( 4 )( 0 ).set( "IQI_name",                   std::string("IQIU_Xtop") );
//     tParameterList( 4 )( 0 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
//     tParameterList( 4 )( 0 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UY ) );
//     tParameterList( 4 )( 0 ).set( "vectorial_field_index",      0 );
//     tParameterList( 4 )( 0 ).set( "master_dof_dependencies",    std::string("UX,UY") );
//     tParameterList( 4 )( 0 ).set( "mesh_set_names",             std::string("iside_g_0_b0_5_b1_1") );
//     
//     tParameterList( 4 )( 1 ) = prm::create_IQI_parameter_list();
//     tParameterList( 4 )( 1 ).set( "IQI_name",                   std::string("IQIU_Xbot") );
//     tParameterList( 4 )( 1 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
//     tParameterList( 4 )( 1 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UY ) );
//     tParameterList( 4 )( 1 ).set( "vectorial_field_index",      1 );
//     tParameterList( 4 )( 1 ).set( "master_dof_dependencies",    std::string("UX,UY") );
//     tParameterList( 4 )( 1 ).set( "mesh_set_names",             std::string("iside_g_0_b0_1_b1_5") );
//     
//     tParameterList( 4 )( 2 ) = prm::create_IQI_parameter_list();
//     tParameterList( 4 )( 2 ).set( "IQI_name",                   std::string("IQIU_Ytop") );
//     tParameterList( 4 )( 2 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
//     tParameterList( 4 )( 2 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UY ) );
//     tParameterList( 4 )( 2 ).set( "vectorial_field_index",      2 );
//     tParameterList( 4 )( 2 ).set( "master_dof_dependencies",    std::string("UX,UY") );
//     tParameterList( 4 )( 2 ).set( "mesh_set_names",             std::string("iside_g_0_b0_5_b1_1") );
//     
//     tParameterList( 4 )( 3 ) = prm::create_IQI_parameter_list();
//     tParameterList( 4 )( 3 ).set( "IQI_name",                   std::string("IQIU_Ybot") );
//     tParameterList( 4 )( 3 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
//     tParameterList( 4 )( 3 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UY ) );
//     tParameterList( 4 )( 3 ).set( "vectorial_field_index",      3 );
//     tParameterList( 4 )( 3 ).set( "master_dof_dependencies",    std::string("UX,UY") );
//     tParameterList( 4 )( 3 ).set( "mesh_set_names",             std::string("iside_g_0_b0_1_b1_5") );


}

void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{
    tParameterlist.resize( 7 );
    for( uint Ik = 0; Ik < 7; Ik ++)
    {
        tParameterlist( Ik ).resize(1);
    }

    tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
    tParameterlist( 0 )(0).set("AZ_diagnostics"    , AZ_none  );
    tParameterlist( 0 )(0).set("AZ_output"         , AZ_none  );
    tParameterlist( 0 )(0).set("AZ_max_iter"       , 10000    );
    tParameterlist( 0 )(0).set("AZ_solver"         , AZ_gmres );
    tParameterlist( 0 )(0).set("AZ_subdomain_solve", AZ_ilu   );
    tParameterlist( 0 )(0).set("AZ_graph_fill"     , 10       );
    tParameterlist( 0 )(0).set("Use_ML_Prec"       ,  false    );

    tParameterlist( 1 )(0) = moris::prm::create_linear_solver_parameter_list();
    tParameterlist( 2 )(0) = moris::prm::create_nonlinear_algorithm_parameter_list();
    tParameterlist( 3 )(0) = moris::prm::create_nonlinear_solver_parameter_list();
    tParameterlist( 3 )(0).set("NLA_DofTypes"      , std::string("UX,UY") );

    tParameterlist( 4 )(0) = moris::prm::create_time_solver_algorithm_parameter_list();
    tParameterlist( 5 )(0) = moris::prm::create_time_solver_parameter_list();
    tParameterlist( 5 )(0).set("TSA_DofTypes"       , std::string("UX,UY") );
    tParameterlist( 5 )(0).set("TSA_Output_Indices" , std::string("0") ); 
    tParameterlist( 5 )(0).set("TSA_Output_Crteria" , std::string("Output_Criterion") );    
    
    tParameterlist( 6 )(0) = moris::prm::create_solver_warehouse_parameterlist();
    
}

void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{       
    tParameterlist.resize( 1 );
    tParameterlist(0).resize(1);
            
    tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    
}

void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{       
    tParameterlist.resize( 1 );
    tParameterlist(0).resize(1);
            
    tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
    tParameterlist( 0 )( 0 ).set( "Mesh_Type"   , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
    tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "/home/sonne/Documents/aaaaa_thesis/benchmark04/" , "benchmark04.exo" ) );
    tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1,HMR_dummy_n_p3,HMR_dummy_c_p3,HMR_dummy_n_p4,HMR_dummy_c_p4,HMR_dummy_n_p5,HMR_dummy_c_p5,HMR_dummy_n_p7,HMR_dummy_c_p7" ) );
    tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY" ) );
    tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL,NODAL" ) );
    tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "UX,UY" ) );

}

void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{       
    tParameterlist.resize( 1 );
    tParameterlist(0).resize(1);
            
    tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

    tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string( "42, 42"));
    tParameterlist( 0 )( 0 ).set( "domain_dimensions", std::string("10, 10") );
    tParameterlist( 0 )( 0 ).set( "domain_offset", std::string("0.0, 0.0") );
    tParameterlist( 0 )( 0 ).set( "domain_sidesets", std::string("1,2,3,4") );
    
    tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",std::string( "0") );
    tParameterlist( 0 )( 0 ).set( "lagrange_input_meshes", std::string("0") );

    tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::string("1" ));
    tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string("0" ));
    tParameterlist( 0 )( 0 ).set( "bspline_orders", std::string("1" ));
    tParameterlist( 0 )( 0 ).set( "bspline_pattern", std::string("0" ));

    tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

    tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
    tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
    tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );
    tParameterlist( 0 )( 0 ).set( "initial_refinement", 0 );

    tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
    tParameterlist( 0 )( 0 ).set( "severity_level", 2 );    
    
    tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 4 );    
}

void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
{       
    tParameterlist.resize( 1 );
    tParameterlist(0).resize(1);
            
    tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

    tParameterlist( 0 )( 0 ).set( "geometries", std::string( "Lvl_set_0,Lvl_set_1,Lvl_set_2") );    
}



//------------------------------------------------------------------------------
        }

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
