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

#include "cl_Param_List.hpp"
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

            // parameter list to build the fem model
            Vector< Vector< ParameterList > > mParameterList;

            // unpacked fem inputs
            Vector< fem::Set_User_Info > mSetInfo;

            // space dimension
            uint mSpaceDim;

            Vector< std::shared_ptr< fem::Field > > mFields;
            Vector< moris::sint >                   mFieldTypes;
            Vector< std::shared_ptr< fem::IQI > >   mIQIs;

            //! requested IQI Names
            Vector< std::string > mRequestedIQINames;

            // A map for mIQIs
            moris::map< std::string, moris_index > mIQINameToIndexMap;

            // flag to skip GEN procedures
            bool mFEMOnly = false;


          public:
            //! Gauss point information. Only used for output
            uint mBulkGaussPoints                 = 0;
            uint mSideSetsGaussPoints             = 0;
            uint mDoubleSidedSideSetsGaussPoints  = 0;
            uint mNonconformalSideSetsGaussPoints = 0;

            /**
             * @brief constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    Vector< fem::Set_User_Info >        &aSetInfo );

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
                    Vector< fem::Set_User_Info >        &aSetInfo,
                    MSI::Design_Variable_Interface      *aDesignVariableInterface );

            /**
             * @brief constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager >     aMeshManager,
                    const moris_index                       &aMeshPairIndex,
                    const Vector< Vector< ParameterList > > &aParameterList,
                    std::shared_ptr< Library_IO >            aLibrary );

            /**
             * @brief constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aLibrary       a file path for property functions
             * @param[ in ] aDesignVariableInterface a design variable interface pointer
             */
            FEM_Model(
                    std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                    const moris_index                   &aMeshPairIndex,
                    Vector< Vector< ParameterList > >    aParameterList,
                    MSI::Design_Variable_Interface      *aDesignVariableInterface );

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
                mBulkGaussPoints                 = 0;
                mSideSetsGaussPoints             = 0;
                mDoubleSidedSideSetsGaussPoints  = 0;
                mNonconformalSideSetsGaussPoints = 0;
            };


            /**
             * @brief resets model member variables
             */
            inline void
            report_on_assembly() override
            {
                uint tTotalBulkGaussPoints                 = sum_all( mBulkGaussPoints );
                uint tTotalSideSetsGaussPoints             = sum_all( mSideSetsGaussPoints );
                uint tTotalDoubleSidedSideSetsGaussPoints  = sum_all( mDoubleSidedSideSetsGaussPoints );
                uint tTotalNonconformalSideSetsGaussPoints = sum_all( mNonconformalSideSetsGaussPoints );

                if ( tTotalBulkGaussPoints + tTotalSideSetsGaussPoints + tTotalDoubleSidedSideSetsGaussPoints + tTotalNonconformalSideSetsGaussPoints > 0 )
                {
                    MORIS_LOG_SPEC( "Number of Bulk Gauss Points", tTotalBulkGaussPoints );
                    MORIS_LOG_SPEC( "Number of SideSet Gauss Points", tTotalSideSetsGaussPoints );
                    MORIS_LOG_SPEC( "Number of DoubleSidedSideset Gauss Points", tTotalDoubleSidedSideSetsGaussPoints );
                    MORIS_LOG_SPEC( "Number of NonconformalSideSet Gauss Points", tTotalNonconformalSideSetsGaussPoints );
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
                    const Vector< gen::PDV_Type > &aPdvTypes,
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
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
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
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
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
                    const Vector< gen::PDV_Type > &aRequestedPdvTypes,
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
                    moris_index        aNodeIndex,
                    enum gen::PDV_Type aPdvType,
                    moris_index        aXYZPdvAssemblyIndex ) override;

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
                    mtk::Interpolation_Mesh      *aIPMesh,
                    mtk::Integration_Mesh        *aIGMesh,
                    Vector< fem::Set_User_Info > &aSetInfo );

            /**
             * @brief Method to update the fem sets that need to be reinitialized in every newton iteration
             */
            void update_equation_sets() override;

            /**
             * @brief set parameter list
             * @param[ in ] aParameterList a list of parameter for the FEM model
             */
            void
            set_parameter_list( const Vector< Vector< ParameterList > > &aParameterList )
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
            set_requested_IQI_names( const Vector< std::string > &aRequestedIQINames ) override
            {
                mRequestedIQINames = aRequestedIQINames;
            }

            /**
             * @brief get requested IQI names
             */
            const Vector< std::string > &
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
             * @brief scale the IQIs according to user input.
             */
            void normalize_IQIs() override;

            /**
             * @brief return field by type
             */
            const std::shared_ptr< fem::Field > &get_field( mtk::Field_Type tFieldType );

            /**
             * @brief return fields
             */
            Vector< std::shared_ptr< mtk::Field > > get_fields() override;


            void populate_fields() override;

            /**
             * @brief get vertex active flags (if relevant) based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aPdvTypes
             * @param[ out ] aIsActiveDv
             */
            void get_vertex_xyz_active_flags( moris_index aVeretxIndex,
                    Matrix< DDSMat >                     &aIsActiveDv,
                    const Vector< enum gen::PDV_Type >   &aPdvTypes );

            /**
             * @brief set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ out ] aIsActiveDv
             */
            void set_vertex_xyz_active_flags( moris_index aVeretxIndex,
                    Vector< Vector< bool > >             &aIsActiveDv );

            /**
             * @brief set vertex active flags (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPvIds list of xyz pdv ids
             */
            void set_vertex_xyz_pdv_ids( moris_index aVeretxIndex,
                    Vector< Matrix< DDSMat > >      &aXYZPvIds );

            /**
             * @brief get vertex xyz pdv ids (if relevant)  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZPdvIds
             * @param[ in ] aIsActiveDv
             */
            void get_vertex_xyz_pdv_ids( moris_index    aVeretxIndex,
                    Matrix< DDSMat >                   &aXYZPdvIds,
                    const Vector< enum gen::PDV_Type > &aPdvTypes );

            /**
             * @brief get x/y/z pdv local cluster assembly indices  based on the index of the vertex
             * @param[ in ] aVeretxIndex
             * @param[ in ] aXYZLocalAssemblyIndices matrix to fill with local cluster assembly indices
             * @param[ in ] aPdvType               list of enums for requested x/y/z pdv
             */
            void get_local_xyz_pdv_assembly_indices( moris_index aVeretxIndex,
                    Matrix< DDSMat >                            &aXYZLocalAssemblyIndices,
                    const Vector< enum gen::PDV_Type >          &aPdvTypes );

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
        };
    }    // namespace fem
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_ */
