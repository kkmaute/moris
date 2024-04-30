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

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"

#include "fn_PRM_FEM_Parameters.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "GEN_Data_Types.hpp"

#include "cl_MSI_Equation_Model.hpp"
#include "cl_FEM_Phase_User_Info.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{

    //------------------------------------------------------------------------------
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
        //------------------------------------------------------------------------------

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
            Vector< fem::Node_Base * > mIPNodes;

            // list of IG node pointers
            Vector< fem::Node_Base * > mIGNodes;

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
            Vector< moris::real > mQi;

            // parameter list to build the fem model
            Vector< Vector< Parameter_List > > mParameterList;

            // unpacked fem inputs
            Vector< fem::Set_User_Info > mSetInfo;

            // unpacked phase inputs
            Vector< fem::Phase_User_Info > mPhaseInfo;
            std::map< std::string, uint >       mPhaseMap;

            // space dimension
            uint mSpaceDim;

            // fixme remove ?
            Vector< std::shared_ptr< fem::Property > >                mProperties;
            Vector< std::shared_ptr< fem::Field > >                   mFields;
            Vector< std::shared_ptr< fem::Material_Model > >          mMMs;
            Vector< std::shared_ptr< fem::Constitutive_Model > >      mCMs;
            Vector< std::shared_ptr< fem::Stabilization_Parameter > > mSPs;
            Vector< std::shared_ptr< fem::IWG > >                     mIWGs;
            Vector< std::shared_ptr< fem::IQI > >                     mIQIs;

            Vector< moris::sint > mFieldTypeMap;

            //! requested IQI Names
            Vector< std::string > mRequestedIQINames;

            // A map for mIQIs
            moris::map< std::string, moris_index > mIQINameToIndexMap;

            // flag to skip GEN procedures
            bool mFEMOnly = false;

            //------------------------------------------------------------------------------

          public:
            //! Gauss point information. Only used for output
            uint mBulkGaussPoints                = 0;
            uint mSideSetsGaussPoints            = 0;
            uint mDoubleSidedSideSetsGaussPoints = 0;

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    Vector< fem::Set_User_Info >   &aSetInfo );

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             * @param[ in ] aDesignVariableInterface a design variable interface pointer
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    Vector< fem::Set_User_Info >   &aSetInfo,
                    MSI::Design_Variable_Interface      *aDesignVariableInterface );

            //------------------------------------------------------------------------------
            /**
             * constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager >        aMeshManager,
                    const moris_index                          &aMeshPairIndex,
                    Vector< Vector< Parameter_List > > aParameterList,
                    std::shared_ptr< Library_IO >               aLibrary );

            //------------------------------------------------------------------------------
            /**
             * constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             * @param[ in ] aDesignVariableInterface a design variable interface pointer
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager >        aMeshManager,
                    const moris_index                          &aMeshPairIndex,
                    Vector< Vector< Parameter_List > > aParameterList,
                    MSI::Design_Variable_Interface             *aDesignVariableInterface );

            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            FEM_Model(){};

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~FEM_Model();

            //------------------------------------------------------------------------------

            void free_memory();

            //------------------------------------------------------------------------------

            /**
             * initialize the FEM model from parameter lists + create the interpolation nodes & FEM sets
             * @param[ in ] aLibrary       a file path for property functions
             */
            void initialize_from_inputfile( std::shared_ptr< Library_IO > aLibrary );

            //------------------------------------------------------------------------------

            /**
             * initialize the FEM model from parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             */
            void initialize( std::shared_ptr< Library_IO > aLibrary );

            //------------------------------------------------------------------------------

            /**
             * resets model member variables
             */
            void
            reset()
            {
                mBulkGaussPoints                = 0;
                mSideSetsGaussPoints            = 0;
                mDoubleSidedSideSetsGaussPoints = 0;
            };

            //------------------------------------------------------------------------------

            /**
             * resets model member variables
             */
            inline void
            report_on_assembly()
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

            //------------------------------------------------------------------------------
            /**
             * create interpolation nodes
             * @param[ in ] aIPMesh interpolation mesh pointer
             */
            void create_interpolation_nodes( mtk::Interpolation_Mesh *aIPMesh );

            //------------------------------------------------------------------------------
            /**
             * create integration nodes
             * @param[ in ] aIGMesh integration mesh pointer
             */
            void create_integration_nodes( mtk::Integration_Mesh *aIGMesh );

            //------------------------------------------------------------------------------
            /**
             * get integration xyz active flags
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_active_flags(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const Vector< gen::PDV_Type > &aPdvTypes,
                    Matrix< DDSMat >              &aIsActiveDv );

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aXYZPdvIds   matrix to fill with ids for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_ids(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aXYZPdvIds );

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             * @param[ in ] aXYZPdvIds   matrix to fill with ids for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_active_flags_and_ids(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aIsActiveDv,
                    Matrix< DDSMat >              &aXYZPdvIds );

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            void get_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >      &aNodeIndices,
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
                    Matrix< DDSMat >              &aXYZPdvAssemblyIndices );

            //------------------------------------------------------------------------------
            
            Matrix< DDSMat >
            get_XYZ_local_pdv_assembly_map() const override;

            //------------------------------------------------------------------------------
            /**
             * reset integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices list of node indices to reset
             */
            void reset_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat > &aNodeIndices );

            //------------------------------------------------------------------------------
            /**
             * set integration xyz pdv local cluster assembly index
             * @param[ in ] aNodeIndex           node index
             * @param[ in ] aPdvType             enum for pdv type
             * @param[ in ] aXYZPdvAssemblyIndex assembly index for pdv type to set
             */
            void set_integration_xyz_pdv_assembly_index(
                    moris_index   aNodeIndex,
                    enum gen::PDV_Type aPdvType,
                    moris_index   aXYZPdvAssemblyIndex );

            //------------------------------------------------------------------------------
            /**
             * create fem sets
             * @param[ in ] aIPMesh interpolation mesh pointer
             * @param[ in ] aIGMesh integration mesh pointer
             */
            void create_fem_sets(
                    mtk::Interpolation_Mesh *aIPMesh,
                    mtk::Integration_Mesh   *aIGMesh );

            //------------------------------------------------------------------------------
            /**
             * create fem sets
             * @param[ in ] aIPMesh  interpolation mesh pointer
             * @param[ in ] aIGMesh  integration mesh pointer
             * @param[ in ] aSetInfo cell of set user info
             */
            void create_fem_sets(
                    mtk::Interpolation_Mesh           *aIPMesh,
                    mtk::Integration_Mesh             *aIGMesh,
                    Vector< fem::Set_User_Info > &aSetInfo );

            //------------------------------------------------------------------------------
            /**
             * set parameter list
             * @param[ in ] aParameterList a list of parameter for the FEM model
             */
            void
            set_parameter_list( Vector< Vector< Parameter_List > > aParameterList )
            {
                mParameterList = aParameterList;
            }

            //------------------------------------------------------------------------------
            /**
             * set space dimension ( only for UT)
             * @param[ in ] aSpaceDim int for space dimension
             */
            void
            set_space_dim( uint aSpaceDim )
            {
                mSpaceDim = aSpaceDim;
            }

            //------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            Vector< MSI::Equation_Set * > &
            get_equation_sets()
            {
                return mFemSets;
            }

            //------------------------------------------------------------------------------
            /**
             * get equation objects
             */
            Vector< MSI::Equation_Object * > &
            get_equation_objects()
            {
                return mFemClusters;
            }

            //------------------------------------------------------------------------------
            /**
             * MTK set to fem set index map
             */
            map< std::tuple< moris_index, bool, bool >, moris_index > &
            get_mesh_set_to_fem_set_index_map()
            {
                return mMeshSetToFemSetMap;
            }

            //------------------------------------------------------------------------------

            /**
             * set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            void
            set_requested_IQI_names( const Vector< std::string > &aRequestedIQINames )
            {
                mRequestedIQINames = aRequestedIQINames;
            }

            //------------------------------------------------------------------------------
            /**
             * get requested IQI names
             */
            const Vector< std::string > &
            get_requested_IQI_names()
            {
                return mRequestedIQINames;
            }

            //------------------------------------------------------------------------------
            /**
             * build a map for the mIQIs, fills in the mIQINameToIndexMap values
             */
            void create_IQI_map();

            //------------------------------------------------------------------------------
            /**
             * Initialize the IQI with the correct size
             */
            void initialize_IQIs();

            //------------------------------------------------------------------------------
            /**
             * finalize the fem sets
             */
            void finalize_equation_sets(
                    MSI::Model_Solver_Interface *aModelSolverInterface );

            //------------------------------------------------------------------------------
            /**
             * create a list of property pointers
             * @param[ in ] aProperties    a list of property pointers to fill
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             * @param[ in ] aLibrary       a file path for property functions
             */
            void create_properties(
                    std::map< std::string, uint >              &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >        &aDvTypeMap,
                    moris::map< std::string, mtk::Field_Type > &aFieldTypeMap,
                    std::shared_ptr< Library_IO >               aLibrary );

            //------------------------------------------------------------------------------
            /**
             * create a list of field pointers
             * * @param[ in ] aFieldNameToIndexMap  Map which maps the field name to an index
             */
            void create_fields(
                    std::map< std::string, uint > &aFieldMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of material model pointers
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_material_models(
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of constitutive model pointers
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_constitutive_models(
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of stabilization parameter pointers
             * @param[ in ] aSPMap         a map from SP name to index in mSPs to fill
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_stabilization_parameters(
                    std::map< std::string, uint >            &aSPMap,
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of IWG pointers
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             */
            void create_IWGs(
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aSPMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create an IQI
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             */
            void create_IQIs(
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aSPMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create fem set info
             * @param[ in ] aWithPhase FIXME remove just there to overload
             */
            void create_fem_set_info( bool aWithPhase );

            //------------------------------------------------------------------------------
            /**
             * create phase info
             */
            void create_phases();

            //------------------------------------------------------------------------------
            /**
             * scale the IQIs according to user input.
             */
            void normalize_IQIs();

            //------------------------------------------------------------------------------
            /**
             * get mesh set name from input
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
                    std::string                 aLeaderPhaseName,
                    std::string                 aFollowerPhaseName,
                    std::string                 aNeighborPhasesString,
                    std::string                 aSideOrdinalsString,
                    bool                        aIsGhost,
                    Vector< std::string > &aMeshSetNames );

            //------------------------------------------------------------------------------
            // FIXME Old version of FEM inputs to be removed
            //------------------------------------------------------------------------------
            /**
             * create a list of material model pointers
             * @param[ in ] aMMMap         a map from MM name to index in aMMs
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_material_models(
                    std::map< std::string, uint >            &aMMMap,
                    std::map< std::string, uint >            &aPropertyMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of constitutive model pointers
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aPropertyMap   a map from property name to property index
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_constitutive_models(
                    std::map< std::string, uint >            &aCMMap,
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aMMMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of stabilization parameter pointers
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             */
            void create_stabilization_parameters(
                    std::map< std::string, uint >            &aSPMap,
                    std::map< std::string, uint >            &aPropertyMap,
                    std::map< std::string, uint >            &aCMMap,
                    moris::map< std::string, MSI::Dof_Type > &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >      &aDvTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create a list of IWG pointers
             * @param[ in ] aPropertyMap   a map from property name to property
             *                             index in aProperties
             * @param[ in ] aCMMap         a map from CM name to CM index in aCMs
             * @param[ in ] aSPMap         a map from SP name to SP index in aSPs
             * @param[ in ] aFieldMap      a map from SP name to index in aFields
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             * @param[ in ] aFieldTypeMap  a map from std::string to Field_Type
             */
            void create_IWGs(
                    std::map< std::string, uint >              &aIWGMap,
                    std::map< std::string, uint >              &aPropertyMap,
                    std::map< std::string, uint >              &aMMMap,
                    std::map< std::string, uint >              &aCMMap,
                    std::map< std::string, uint >              &aSPMap,
                    std::map< std::string, uint >              &aFieldMap,
                    moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >        &aDvTypeMap,
                    moris::map< std::string, mtk::Field_Type > &aFieldTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create an IQI
             * @param[ in ] aPropertyMap   a map from property name to index in mProperties
             * @param[ in ] aCMMap         a map from CM name to index in aCMs
             * @param[ in ] aSPMap         a map from SP name to index in aSPs
             * @param[ in ] aFieldMap      a map from SP name to index in aFields
             * @param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * @param[ in ] aDvTypeMap     a map from std::string to gen::PDV_Type
             * @param[ in ] aFieldTypeMap  a map from std::string to Field_Type
             */
            void create_IQIs(
                    std::map< std::string, uint >              &aIQIMap,
                    std::map< std::string, uint >              &aPropertyMap,
                    std::map< std::string, uint >              &aCMMap,
                    std::map< std::string, uint >              &aSPMap,
                    std::map< std::string, uint >              &aFieldMap,
                    moris::map< std::string, MSI::Dof_Type >   &aMSIDofTypeMap,
                    moris::map< std::string, gen::PDV_Type >        &aDvTypeMap,
                    moris::map< std::string, mtk::Field_Type > &aFieldTypeMap );

            //------------------------------------------------------------------------------
            /**
             * create fem set info
             */
            void create_fem_set_info();

            /**
             * Updates the stored fields.
             */
            void update_fields() override;

            /**
             * return field by type
             */
            const std::shared_ptr< fem::Field > &get_field( mtk::Field_Type tFieldType );

            //------------------------------------------------------------------------------
            /**
             * return fields
             */
            Vector< std::shared_ptr< mtk::Field > > get_fields();

            //------------------------------------------------------------------------------

            void populate_fields();

            //------------------------------------------------------------------------------
            /**
             * get vertex active flags (if relevant) based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aPdvTypes
             * @param[ out ] aIsActiveDv
             */
            void get_vertex_xyz_active_flags( moris_index aVeretxIndex,
                    Matrix< DDSMat >                     &aIsActiveDv,
                    const Vector< enum gen::PDV_Type >   &aPdvTypes );

            //------------------------------------------------------------------------------
            /**
             * set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aIsActiveDv
             */
            void set_vertex_xyz_active_flags(
                    moris_index                 aVertexIndex,
                    const Vector< Vector< bool > >& aIsActiveDv );

            //------------------------------------------------------------------------------
            /**
             * set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPvIds list of xyz pdv ids
             */
            void set_vertex_xyz_pdv_ids( moris_index aVeretxIndex,
                    Vector< Matrix< DDSMat > > &aXYZPvIds );

            //------------------------------------------------------------------------------
            /**
             * get vertex xyz pdv ids (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPdvIds
             * @param[ in ] aIsActiveDv
             */
            void get_vertex_xyz_pdv_ids( moris_index    aVeretxIndex,
                    Matrix< DDSMat >                   &aXYZPdvIds,
                    const Vector< enum gen::PDV_Type > &aPdvTypes );

            //------------------------------------------------------------------------------
            /**
             * get x/y/z pdv local cluster assembly indices  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZLocalAssemblyIndices matrix to fill with local cluster assembly indices
             * @param[ in ] aPdvType               list of enums for requested x/y/z pdv
             */
            void get_local_xyz_pdv_assembly_indices( moris_index aVeretxIndex,
                    Matrix< DDSMat >                            &aXYZLocalAssemblyIndices,
                    const Vector< enum gen::PDV_Type >          &aPdvTypes );


            //------------------------------------------------------------------------------

            /**
             * @brief Set the dof type to Bspline mesh index map
             *
             * @param aDofTypeToBsplineMeshIndex
             */
            void
            set_dof_type_to_Bspline_mesh_index( std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex );

            //------------------------------------------------------------------------------

            /**
             * @brief set flag whether to use new ghost sets
             *
             * @param aUseNewGhostSets
             */
            void
            set_use_new_ghost_sets( bool aUseNewGhostSets );

            //------------------------------------------------------------------------------

            /**
             * @brief check the ghost set names and conrrect them to automatically include
             * the correct B-spline mesh index when the new ghost sets are used
             *
             * @param aMeshSetName
             * @param aDofType
             */
            void
            check_and_set_ghost_set_names( std::string &aMeshSetName, enum MSI::Dof_Type aDofType );

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    }    // namespace fem
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_ */
