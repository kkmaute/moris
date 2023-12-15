/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Model.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"

#include "fn_PRM_FEM_Parameters.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

#include "cl_MSI_Equation_Model.hpp"
#include "cl_FEM_Phase_User_Info.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{

    namespace mtk
    {
        class Mesh_Manager;
        class Interpolation_Mesh;
        class Integration_Mesh;
        class Field;
        enum class Field_Type;
    }    // namespace mtk

    namespace fem
    {
        class IWG;
        class Node_Base;
        class Set;
        class Field_Interpolator;
        class Property;
        class Material_Model;
        class Constitutive_Model;
        class Stabilization_Parameter;
        class IWG;
        class IQI;
        class Field;
    }    // namespace fem

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
    }    // namespace MSI

    namespace fem
    {

        class FEM_Model : public MSI::Equation_Model
        {
            // pointer to reference mesh
            std::shared_ptr< mtk::Mesh_Manager > mMeshManager = nullptr;
            moris_index                          mMeshPairIndex;

            // information about B-spline mesh indices (only for naming purposes, this information is otherwise entity of MSI)
            std::unordered_map< MSI::Dof_Type, moris_index > mDofTypeToBsplineMeshIndex;

            // flag whether to use the new ghost sets or not
            bool mUseNewGhostSets = false;

            // list of IP node pointers
            moris::Cell< fem::Node_Base * > mIPNodes;

            // list of IG node pointers
            moris::Cell< fem::Node_Base * > mIGNodes;

            // storage for xyz pdv active flags
            // row number i : underlying mtk vertex index of the fem IG node
            // column number : showing activation based on type
            Matrix< DDSMat > mIsActiveXYZ;

            // storage for xyz pdv ids
            // row number i : underlying mtk vertex index of the fem IG node
            // column number : showing activation based on type
            Matrix< DDSMat > mXYZPdvIds;

            // storage for xyz pdv ids
            // row number i : underlying mtk vertex index of the fem IG node
            // column number : showing activation based on type
            Matrix< DDSMat > mXYZLocalAssemblyIndices;

            // list of QI values
            moris::Cell< moris::real > mQi;

            // parameter list to build the fem model
            moris::Cell< moris::Cell< ParameterList > > mParameterList;

            // unpacked fem inputs
            moris::Cell< fem::Set_User_Info > mSetInfo;

            // unpacked phase inputs
            moris::Cell< fem::Phase_User_Info > mPhaseInfo;
            std::map< std::string, uint >       mPhaseMap;

            // space dimension
            uint mSpaceDim;

            // fixme remove ?
            moris::Cell< std::shared_ptr< fem::Property > >                mProperties;
            moris::Cell< std::shared_ptr< fem::Field > >                   mFields;
            moris::Cell< std::shared_ptr< fem::Material_Model > >          mMMs;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      mCMs;
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mSPs;
            moris::Cell< std::shared_ptr< fem::IWG > >                     mIWGs;
            moris::Cell< std::shared_ptr< fem::IQI > >                     mIQIs;

            moris::Cell< moris::sint > mFieldTypeMap;

            //! requested IQI Names
            moris::Cell< std::string > mRequestedIQINames;

            // A map for mIQIs
            moris::map< std::string, moris_index > mIQINameToIndexMap;

            // flag to skip GEN procedures
            bool mFEMOnly = false;


          public:
            //! Gauss point information. Only used for output
            uint mBulkGaussPoints                = 0;
            uint mSideSetsGaussPoints            = 0;
            uint mDoubleSidedSideSetsGaussPoints = 0;

            /**
             * @brief constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    moris::Cell< fem::Set_User_Info >   &aSetInfo );

            /**
             * @brief constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             * @param[ in ] aDesignVariableInterface a design variable interface pointer
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    moris::Cell< fem::Set_User_Info >   &aSetInfo,
                    MSI::Design_Variable_Interface      *aDesignVariableInterface );

            /**
             * @brief constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager >               aMeshManager,
                    const moris_index                                 &aMeshPairIndex,
                    const moris::Cell< moris::Cell< ParameterList > > &aParameterList,
                    std::shared_ptr< Library_IO >                      aLibrary );

            /**
             * @brief constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             * @param[ in ] aDesignVariableInterface a design variable interface pointer
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager >        aMeshManager,
                    const moris_index                          &aMeshPairIndex,
                    moris::Cell< moris::Cell< ParameterList > > aParameterList,
                    MSI::Design_Variable_Interface             *aDesignVariableInterface );

            /**
             * @brief trivial constructor
             */
            FEM_Model() = default;

            ~FEM_Model();

            void free_memory() final;


            /**
             * @brief initialize the FEM model from parameter lists + create the interpolation nodes & FEM sets
             * @param[ in ] aLibrary       a file path for property functions
             */
            void initialize_from_inputfile( std::shared_ptr< Library_IO > aLibrary ) override;


            /**
             * @brief initialize the FEM model from parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             */
            void initialize( const std::shared_ptr< Library_IO > &aLibrary );


            /**
             * @brief resets model member variables
             */
            void
            reset() override
            {
                mBulkGaussPoints                = 0;
                mSideSetsGaussPoints            = 0;
                mDoubleSidedSideSetsGaussPoints = 0;
            };


            /**
             * @brief resets model member variables
             */
            inline void
            report_on_assembly() override
            {
                uint tTotalBulkGaussPoints                = sum_all( mBulkGaussPoints );
                uint tTotalSideSetsGaussPoints            = sum_all( mSideSetsGaussPoints );
                uint tTotalDoubleSidedSideSetsGaussPoints = sum_all( mDoubleSidedSideSetsGaussPoints );

                if ( tTotalBulkGaussPoints + tTotalSideSetsGaussPoints + tTotalDoubleSidedSideSetsGaussPoints > 0 )
                {
                    MORIS_LOG_SPEC( "Number of Bulk Gauss Points", tTotalBulkGaussPoints );
                    MORIS_LOG_SPEC( "Number of SideSet Gauss Points", tTotalSideSetsGaussPoints );
                    MORIS_LOG_SPEC( "Number of DoubleSidedSideset Gauss Points", tTotalDoubleSidedSideSetsGaussPoints );
                }
            };

            /**
             * @brief create interpolation nodes
             * @param[ in ] aIPMesh interpolation mesh pointer
             */
            void create_interpolation_nodes( mtk::Interpolation_Mesh *aIPMesh );

            /**
             * @brief create integration nodes
             * @param[ in ] aIGMesh integration mesh pointer
             */
            void create_integration_nodes( mtk::Integration_Mesh *aIGMesh );

            /**
             * @brief get integration xyz active flags
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_active_flags(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const moris::Cell< PDV_Type > &aPdvTypes,
                    Matrix< DDSMat >              &aIsActiveDv ) override;

            /**
             * @brief get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aXYZPdvIds   matrix to fill with ids for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_ids(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const moris::Cell< PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aXYZPdvIds ) override;

            /**
             * @brief get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             * @param[ in ] aXYZPdvIds   matrix to fill with ids for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_active_flags_and_ids(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const moris::Cell< PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aIsActiveDv,
                    Matrix< DDSMat >              &aXYZPdvIds ) override;

            /**
             * @brief get integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const moris::Cell< PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aXYZPdvAssemblyIndices ) override;


            Matrix< DDSMat >
            get_XYZ_local_pdv_assembly_map() const override;

            /**
             * @brief reset integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices list of node indices to reset
             */
            void reset_integration_xyz_pdv_assembly_indices( const Matrix< IndexMat > &aNodeIndices ) override;

            /**
             * @brief set integration xyz pdv local cluster assembly index
             * @param[ in ] aNodeIndex           node index
             * @param[ in ] aPdvType             enum for pdv type
             * @param[ in ] aXYZPdvAssemblyIndex assembly index for pdv type to set
             */
            void set_integration_xyz_pdv_assembly_index(
                    moris_index   aNodeIndex,
                    enum PDV_Type aPdvType,
                    moris_index   aXYZPdvAssemblyIndex ) override;

            /**
             * @brief create fem sets
             * @param[ in ] aIPMesh interpolation mesh pointer
             * @param[ in ] aIGMesh integration mesh pointer
             */
            void create_fem_sets(
                    mtk::Interpolation_Mesh *aIPMesh,
                    mtk::Integration_Mesh   *aIGMesh );

            /**
             * @brief create fem sets
             * @param[ in ] aIPMesh  interpolation mesh pointer
             * @param[ in ] aIGMesh  integration mesh pointer
             * @param[ in ] aSetInfo cell of set user info
             */
            void create_fem_sets(
                    mtk::Interpolation_Mesh           *aIPMesh,
                    mtk::Integration_Mesh             *aIGMesh,
                    moris::Cell< fem::Set_User_Info > &aSetInfo );

            /**
             * @brief set parameter list
             * @param[ in ] aParameterList a list of parameter for the FEM model
             */
            void
            set_parameter_list( const moris::Cell< moris::Cell< ParameterList > > &aParameterList )
            {
                mParameterList = aParameterList;
            }

            /**
             * @brief set space dimension ( only for UT)
             * @param[ in ] aSpaceDim int for space dimension
             */
            void
            set_space_dim( uint aSpaceDim )
            {
                mSpaceDim = aSpaceDim;
            }


            /**
             * @brief set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            void
            set_requested_IQI_names( const moris::Cell< std::string > &aRequestedIQINames ) override
            {
                mRequestedIQINames = aRequestedIQINames;
            }

            /**
             * @brief get requested IQI names
             */
            const moris::Cell< std::string > &
            get_requested_IQI_names() override
            {
                return mRequestedIQINames;
            }

            /**
             * @brief build a map for the mIQIs, fills in the mIQINameToIndexMap values
             */
            void create_IQI_map() override;

            /**
             * @brief Initialize the IQI with the correct size
             */
            void initialize_IQIs() override;

            /**
             * @brief finalize the fem sets
             */
            void finalize_equation_sets( MSI::Model_Solver_Interface *aModelSolverInterface ) override;

            /**
             * @brief create a list of property pointers
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             * @param[ in ] aLibrary       a file path for property functions
             */
            std::map< std::string, uint > create_properties(moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                    moris::map< std::string, PDV_Type >        &aDvTypeMap,
                    moris::map< std::string, mtk::Field_Type > &aFieldTypeMap,
                    std::shared_ptr< Library_IO >               aLibrary );

            /**
             * @brief create a list of field pointers
             */
            void create_fields( );

            /**
             * @brief create a list of material model pointers
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            void create_material_models(
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, PDV_Type >      &aDvTypeMap );

            /**
             * @brief create a list of constitutive model pointers
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            void create_constitutive_models(
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, PDV_Type >      &aDvTypeMap );

            /**
             * @brief create a list of stabilization parameter pointers
             * @param[ in ] aSPMap         a map from SP name to index in mSPs to fill
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            std::map< std::string, uint > create_stabilization_parameters( std::map< std::string, uint > &aPropertyMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap );

            /**
             * @brief create a list of IWG pointers
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             */
            void create_IWGs(
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aSPMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap );

            /**
             * @brief create an IQI
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             */
            void create_IQIs(
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aSPMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap );

            /**
             * @brief create fem set info
             */
            void create_fem_set_info();

            /**
             * @brief create phase info
             */
            void create_phases();

            /**
             * @brief scale the IQIs according to user input.
             */
            void normalize_IQIs() override;

            /**
             * @brief get mesh set name from input
             * @param[ in ] aBulkType             enum for bulk type
             *                                    (bulk, single sideset, double sideset, ...)
             * @param[ in ] aLeaderPhaseName      name for leader phase
             * @param[ in ] aFollowerPhaseName       name for follower phase
             * @param[ in ] aNeighborPhasesString string with neighboring phases for single sideset
             * @param[ in ] aSideOrdinalsString   string with side ordinals for single sideset
             * @param[ in ] aIsGhost              bool true if ghost IWG
             * @param[ in ] aMeshSetNames         cell of mesh set names to fill
             */
            void get_mesh_set_names(
                    fem::Element_Type           aBulkType,
                    const std::string          &aLeaderPhaseName,
                    const std::string          &aFollowerPhaseName,
                    const std::string          &aNeighborPhasesString,
                    const std::string          &aSideOrdinalsString,
                    bool                        aIsGhost,
                    moris::Cell< std::string > &aMeshSetNames );

            // FIXME Old version of FEM inputs to be removed
            /**
             * @brief create a list of material model pointers
             * @param[ in ] aMMMap         a map from MM name to index in aMMs
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            std::map< std::string, uint > create_material_models_without_phase( std::map< std::string, uint > &aPropertyMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap );

            /**
             * @brief create a list of constitutive model pointers
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            std::map< std::string, uint > create_constitutive_models_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aMMMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap );

            /**
             * @brief create a list of stabilization parameter pointers
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             */
            std::map< std::string, uint > create_stabilization_parameters_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aCMMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap );

            /**
             * @brief create a list of IWG pointers
             * @param[ in ] aPropertyMap   a map from property name to property
             *                             index in aProperties
             * @param[ in ] aCMMap         a map from CM name to CM index in aCMs
             * @param[ in ] aSPMap         a map from SP name to SP index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             * @param[ in ] aFieldTypeMap  a map from std::string to Field_Type
             */
            void create_IWGs_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aMMMap, std::map< std::string, uint > &aCMMap, std::map< std::string, uint > &aSPMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap, moris::map< std::string, mtk::Field_Type > &aFieldTypeMap );

            /**
             * @brief create an IQI
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to PDV_Type
             * @param[ in ] aFieldTypeMap  a map from std::string to Field_Type
             */
            void create_IQIs_without_phase( std::map< std::string, uint > &aPropertyMap, std::map< std::string, uint > &aCMMap, std::map< std::string, uint > &aSPMap, moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, moris::map< std::string, PDV_Type > &aDvTypeMap, moris::map< std::string, mtk::Field_Type > &aFieldTypeMap );

            /**
             * @brief create fem set info
             */
            void create_fem_set_info_without_phase();

            /**
             * @brief return field by type
             */
            const std::shared_ptr< fem::Field > &get_field( mtk::Field_Type tFieldType );

            /**
             * @brief return fields
             */
            moris::Cell< std::shared_ptr< mtk::Field > > get_fields() override;


            void populate_fields() override;

            /**
             * @brief get vertex active flags (if relevant) based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aPdvTypes
             * @param[ out ] aIsActiveDv
             */
            void get_vertex_xyz_active_flags( moris_index aVeretxIndex,
                    Matrix< DDSMat >                     &aIsActiveDv,
                    const moris::Cell< enum PDV_Type >   &aPdvTypes );

            /**
             * @brief set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aIsActiveDv
             */
            void set_vertex_xyz_active_flags( moris_index aVeretxIndex,
                    moris::Cell< Matrix< DDSMat > >      &aIsActiveDv );

            /**
             * @brief set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPvIds list of xyz pdv ids
             */
            void set_vertex_xyz_pdv_ids( moris_index aVeretxIndex,
                    moris::Cell< Matrix< DDSMat > > &aXYZPvIds );

            /**
             * @brief get vertex xyz pdv ids (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPdvIds
             * @param[ in ] aIsActiveDv
             */
            void get_vertex_xyz_pdv_ids( moris_index    aVeretxIndex,
                    Matrix< DDSMat >                   &aXYZPdvIds,
                    const moris::Cell< enum PDV_Type > &aPdvTypes );

            /**
             * @brief get x/y/z pdv local cluster assembly indices  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZLocalAssemblyIndices matrix to fill with local cluster assembly indices
             * @param[ in ] aPdvType               list of enums for requested x/y/z pdv
             */
            void get_local_xyz_pdv_assembly_indices( moris_index aVeretxIndex,
                    Matrix< DDSMat >                            &aXYZLocalAssemblyIndices,
                    const moris::Cell< enum PDV_Type >          &aPdvTypes );


            /**
             * @brief Set the dof type to Bspline mesh index map
             * @param aDofTypeToBsplineMeshIndex
             */
            void
            set_dof_type_to_Bspline_mesh_index(
                    std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex ) override;


            /**
             * @brief set flag whether to use new ghost sets
             * @param aUseNewGhostSets
             */
            void
            set_use_new_ghost_sets( bool aUseNewGhostSets ) override;


            /**
             * @brief check the ghost set names and conrrect them to automatically include
             * the correct B-spline mesh index when the new ghost sets are used
             * @param aMeshSetName
             * @param aDofType
             */
            void
            check_and_set_ghost_set_names( std::string &aMeshSetName, enum MSI::Dof_Type aDofType );


          private:
            void set_IWG_dof_dependencies( map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower aLeaderFollowerType ) const;
            void set_IWG_dv_dependencies( map< std::string, PDV_Type > &aDvTypeMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
            void set_IWG_field_types( map< std::string, mtk::Field_Type > &aFieldTypeMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const;
            void set_IWG_properties( std::map< std::string, uint > &aPropertyMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
            void set_IWG_material_models( std::map< std::string, uint > &aMMMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
            void set_IWG_constitutive_models( std::map< std::string, uint > &aCMMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType );
            void set_IWG_stabilization_parameters( std::map< std::string, uint > &aSPMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG );
            void set_IWG_residual_dof_type( map< std::string, MSI::Dof_Type > &aMSIDofTypeMap, ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG ) const;
            void create_fem_set_info_from_IWGs( bool aIsAnalyticalSA, FDScheme_Type const &aFDSchemeForSA, real const aFDPerturbation, std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );
            void create_fem_set_info_from_IQIs( bool const aIsAnalyticalSA, FDScheme_Type const &aFDSchemeForSA, real const aFDPerturbation, std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex );


            template< typename T >
            moris::Cell< moris::Cell< T > >
            property_to_cell_of_cell( ParameterList const &aIWGParameter, std::string const &aPropertyName, map< std::string, T > const &aTypeMap ) const
            {
                return string_to_cell_of_cell< T >( aIWGParameter.get< std::string >( aPropertyName ), aTypeMap );
            }
            void print_physics_model( bool aWithPhase );
        };
    }    // namespace fem
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_ */
