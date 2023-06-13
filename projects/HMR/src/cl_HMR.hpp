/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_HPP_
#define SRC_HMR_CL_HMR_HPP_

//#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Element.hpp"
#include "cl_HMR_Field_Param.hpp"
#include "cl_HMR_Parameters.hpp"     //HMR/src
#include "cl_Library_IO.hpp"

namespace moris::mtk
{
    class Mesh_Manager;
}
namespace moris::hmr
{
    class Field;
    class Database;
    class Mesh;
    class Lagrange_Mesh_Base;
    class Interpolation_Mesh_HMR;
    class Integration_Mesh_HMR;
    // -----------------------------------------------------------------------------
    /**
     * \brief the main class of HMR
     */
    class HMR
    {
        public :
            //! object containing user settings
            Parameters * mParameters;

            std::shared_ptr< Database > mDatabase;

            //! flag telling if perform_refinement() has been called
            bool mUpdateRefinementCalled  = false;

            bool mRestartedFromFile = false;

            //! mesh which points to input pattern
            Cell< std::shared_ptr< Mesh > > mMeshes;
            Cell< std::shared_ptr< Mesh > > mInputMeshes;

            //! mesh which points to output pattern
            //Cell< std::shared_ptr< Mesh > > mOutputMeshes;

            //! container with field objects
            Cell< std::shared_ptr< Field > > mFields;

            //! map for Lagrange orders
            Matrix< DDUMat > mLagrangeOrderToInputMeshIndexMap;

            std::shared_ptr< mtk::Mesh_Manager > mMTKPerformer = nullptr;

            /*
             * @brief determines elements (cells) intersected by the level set
             *
             * @param[in] aCells        - elements to be flagged for refinement
             * @param[in] aCandidates   - candidates for refinement
             * @param[in] aVertexValues - vertex values of scalar field
             * @param[in] aLowerBound   - lower bound of LS
             * @param[in] aUpperBound   - upper bound of LS
             */
            void find_cells_intersected_by_levelset(
                    Cell< hmr::Element * >   & aCells,
                    Cell< hmr::Element * >   & aCandidates,
                    const  Matrix< DDRMat >  & aVertexValues,
                    real                       aLowerBound = -0.0001,
                    real                       aUpperBound =  0.0001);

            void find_low_level_cells_intersected_by_levelset(
                    Cell< hmr::Element * >   & aCells,
                    Cell< hmr::Element * >   & aCandidates,
                    const  Matrix< DDRMat >  & aVertexValues,
                    real                       aLowerBound = -0.0001,
                    real                       aUpperBound =  0.0001);

            /*
             * @brief determines volume elements (cells)
             *
             * @param[in] aCells        - elements to be flagged for refinement
             * @param[in] aCandidates   - candidates for refinement
             * @param[in] aVertexValues - vertex values of scalar field
             * @param[in] aUpperBound   - upper bound of LS
             */
            void find_cells_within_levelset(
                    Cell< hmr::Element * >  & aCells,
                    Cell< hmr::Element * >  & aCandidates,
                    const  Matrix< DDRMat > & aVertexValues,
                    uint                      aUpperBound = 0.0 );

            // -----------------------------------------------------------------------------
        public :
            // -----------------------------------------------------------------------------

            /**
             * the  HMR constructor
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR ( Parameters * aParameters ) ;

            // -----------------------------------------------------------------------------

            /**
             * alternative constructor using a ref
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR ( Parameters & aParameters ) ;

            // -----------------------------------------------------------------------------

            /**
             * alternative constructor using a parameter list
             *
             * @param[in] aParameters  ref to container of user defined settings
             */
            HMR (
                    ParameterList                       & aParameterList,
                    std::shared_ptr<moris::Library_IO>    aLibrary = nullptr ) ;

            // -----------------------------------------------------------------------------

            /**
             * alternative constructor which loads a mesh from a h5 file
             */
            HMR( const std::string & aPath );

            // -----------------------------------------------------------------------------

            /**
             * alternative constructor which loads input and output patterns from path
             */
            HMR(
                    const std::string & aInPath,
                    const std::string & aOutPath );

            // -----------------------------------------------------------------------------

            /**
             * default destructor of HMR
             */
            ~HMR ();

            // -----------------------------------------------------------------------------
            /**
             * save the mesh to an exodus file
             */
            void save_to_exodus(
                    uint aMeshIndex,
                    const std::string & aPath,
                    double              aTimeStep = 0.0 );

            // -----------------------------------------------------------------------------

            /**
             * save the mesh to an exodus file
             */
            void save_last_step_to_exodus(
                    uint                aIndex,
                    const std::string & aPath,
                    double              aTimeStep = 0.0 );

            // -----------------------------------------------------------------------------

            /**
             * save the mesh to an hdf5 file
             */
            void save_to_hdf5(
                    const std::string & aPath,
                    uint                aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            /**
             * store the T-Matrices and B-Spline IDs into a file
             */
            void save_coeffs_to_hdf5_file(
                    const std::string & aFilePath,
                    uint aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            /**
             * store the T-Matrices and B-Spline IDs into a file
             */
            void save_mesh_relations_to_hdf5_file (
                    const std::string & aFilePath,
                    uint aLagrangeMeshIndex,
                    uint aBsplineMeshIndex );

            // -----------------------------------------------------------------------------

            /**
             * loads a field from an HDF5 file and creates a smart pointer
             * to it
             */
            std::shared_ptr< Field > load_field_from_hdf5_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                          uint          aLagrangeIndex = 0,
                          uint          aBSpineIndex = 0 );

            // -----------------------------------------------------------------------------

            std::shared_ptr< Field > load_field_from_exo_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                          uint          aLagrangeIndex=0,
                          uint          aBSpineIndex=0 );

            // -----------------------------------------------------------------------------

            std::shared_ptr< Field > load_field_from_file(
                    const std::string & aLabel,
                    const std::string & aFilePath,
                          uint          aLagrangeIndex=0,
                          uint          aBSpineIndex=0 );

            // -----------------------------------------------------------------------------

            // -----------------------------------------------------------------------------

            /**
             * flags active elements on working pattern
             *
             * @param[ in ]   aElements            element pointers that are to be flagged
             * @param[ in ]   aMinRefinementLevel  if the level of the child is less than this value
             *                                     the child is automatically flagged in the next iteration
             */
            void flag_elements_on_working_pattern(
                    Cell< hmr::Element* > & aElements,
                    uint               aMinRefinementLevel = 0 );

            // -----------------------------------------------------------------------------

            /**
             * exposes the parameters pointer
             */
            Parameters * get_parameters()
            {
                return mParameters;
            }

            // -----------------------------------------------------------------------------

            /**
             * runs the refinement scheme
             */
            void perform_refinement_based_on_working_pattern(
                    uint aPattern,
                    bool aResetPattern = false );

            void perform_refinement( uint aPattern );

            // -----------------------------------------------------------------------------

            void put_elements_on_refinment_queue( Cell< hmr::Element* > & aElements);

            // -----------------------------------------------------------------------------

            /**
             * copy output pattern to input pattern
             */
            void update_refinement_pattern( uint aPattern) ;

            // -----------------------------------------------------------------------------

            /**
             * set activation pattern
             */
            void set_activation_pattern( uint aActivationPattern );

            // -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Max Lagrange Order, Outpot pattern
             */
            std::shared_ptr< Mesh > create_mesh();

            // -----------------------------------------------------------------------------

            /**
             * Creates an STK interface object.
             * Default: Output pattern
             */
            std::shared_ptr< Mesh > create_mesh( uint aLagrangeOrder );

            // -----------------------------------------------------------------------------

            std::shared_ptr< Mesh > create_mesh( uint aLagrangeOrder,
                    uint aPattern );

            std::shared_ptr< Mesh > create_mesh(
                    uint aLagrangeOrder,
                    uint aLagrangePattern,
                    uint aBsplinePattern );

            // -----------------------------------------------------------------------------

            Interpolation_Mesh_HMR * create_interpolation_mesh( uint aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            Interpolation_Mesh_HMR * create_interpolation_mesh(
                    const std::string & aName,
                    bool                aTMatricesExist = false );

            // -----------------------------------------------------------------------------

            Interpolation_Mesh_HMR * create_interpolation_mesh(
                    uint aLagrangeOrder,
                    uint aPattern );

            Interpolation_Mesh_HMR * create_interpolation_mesh(
                    uint aOrder,
                    uint aLagrangePattern,
                    uint aBsplinePattern);

            // -----------------------------------------------------------------------------

            Integration_Mesh_HMR * create_integration_mesh(
                    uint aLagrangeOrder,
                    uint aPattern,
                    Interpolation_Mesh_HMR * aInterpolationMesh);

            // -----------------------------------------------------------------------------

            Integration_Mesh_HMR * create_integration_mesh(
                    uint aLagrangeMeshIndex,
                    Interpolation_Mesh_HMR * aInterpolationMesh);

            // -----------------------------------------------------------------------------
            std::shared_ptr< Field > create_field( const std::string & aLabel );

            // -----------------------------------------------------------------------------

            std::shared_ptr< Field > create_field(
                    const std::string & aLabel,
                    uint aLagrangeOrder,
                    uint aBSplineOrder );

            // -----------------------------------------------------------------------------

            // create field from parameter list
            // std::shared_ptr< Field > create_field( const Field_Param & aParameters );

            // -----------------------------------------------------------------------------

            /**
             * grab the pointer to the datavase
             */
            // std::shared_ptr< Database >
            auto get_database() -> decltype ( mDatabase )
            {
                return mDatabase;
            }

            // -----------------------------------------------------------------------------

            /**
             * Flag an element for refinement. Needed for Testing.
             */
            void flag_element( moris_index aElementIndex );

            // -----------------------------------------------------------------------------

            /**
             * for flagging
             */
            void get_candidates_for_refinement(
                    Cell< hmr::Element* > & aCandidates,
                    uint              aLagrangeMeshIndex);

            void get_candidates_for_refinement(
                    Cell< hmr::Element* > & aCandidates,
                    Lagrange_Mesh_Base    * aMesh );

            void get_active_candidates_for_refinement(
                    Cell< hmr::Element* > & aCandidates,
                    uint              aLagrangeMeshIndex);

            void get_active_candidates_for_refinement(
                    Cell< hmr::Element* > & aCandidates,
                    Lagrange_Mesh_Base    * aMesh );

            // -----------------------------------------------------------------------------

            /**
             * flags elements on the surface and inside of a level set
             */
            uint flag_volume_and_surface_elements_on_working_pattern( std::shared_ptr<Field> aScalarField );

            // -----------------------------------------------------------------------------

            /**
             * flags elements on the surface of a level set
             */
            uint flag_surface_elements_on_working_pattern( std::shared_ptr<Field> aScalarField );

            //FIXME this can be deleted the moment GEN is fixed
            uint based_on_field_put_elements_on_queue(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aLagrangeMeshIndex,
                    sint                     aFunctionIndex = -1);

            uint based_on_field_put_elements_on_queue(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aPattern,
                    uint                     aOrder,
                    Refinement_Function      aRefFunction);

            uint based_on_field_flag_elements_on_working_pattern(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aPattern,
                    uint                     aOrder,
                    Refinement_Function      aRefFunction );

            uint based_on_field_put_low_level_elements_on_queue(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aLagrangeMeshIndex,
                    sint                     aFunctionIndex);

            uint based_on_field_put_low_level_elements_on_queue(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aPattern,
                    uint                     aOrder,
                    sint                     aFunctionIndex);

            uint based_on_field_put_low_level_elements_on_queue(
                    const Matrix< DDRMat > & aFieldValues,
                    uint                     aPattern,
                    uint                     aOrder,
                    Refinement_Function      aRefFunction);

            // -----------------------------------------------------------------------------

            /**
             * special function for tutorial
             */
            //            void perform_refinement_and_map_fields( const uint aPattern );

            // -----------------------------------------------------------------------------

            /**
             * funciton for L2 test
             */
            void map_field_to_output(
                    std::shared_ptr< Field > aField,
                    uint                     aMesh_Index,
                    uint                     aBsplineMeshIndex);

            void map_field_to_output_union(
                    std::shared_ptr< Field > aField,
                    uint                     aUnionOrder );

            // -----------------------------------------------------------------------------

            /**
             * calculate T-Matrices, faces and edges
             */
            void finalize();

            // -----------------------------------------------------------------------------

            /**
             * resets HMR. deletes T-Matrices, faces and edges
             */

            void reset_HMR();

            // -----------------------------------------------------------------------------

            /**
             * performs HMR specific task. Right now only finalize. more will be added
             */
            void perform();

            // -----------------------------------------------------------------------------

            void set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

            // -----------------------------------------------------------------------------

            bool get_mesh_name_exists( const std::string & aName ) const;

            // -----------------------------------------------------------------------------

            void load_refinement_from_file();

            // -----------------------------------------------------------------------------

            void output_mesh_refinement_data();

            // -----------------------------------------------------------------------------
            // Debug files
            // -----------------------------------------------------------------------------

            /**
             * Dumps the background mesh into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void save_background_mesh_to_vtk( const std::string & aFilePath );

            // -----------------------------------------------------------------------------

            /**
             * Dumps the Bsplines into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void save_bsplines_to_vtk(
                    const std::string& aFilePath,
                    uint               aLagrangeMeshIndex,
                    uint               aBsplineMeshIndex  );

            // -----------------------------------------------------------------------------

            void calculate_bspline_coordinates(
                    uint aLagrangeMeshIndex,
                    uint aBsplineMeshIndex  );

            // -----------------------------------------------------------------------------

            /**
             * Dumps the faces into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void save_faces_to_vtk(
                    const std::string & aFilePath,
                    uint aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            /**
             * Dumps the edges into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void save_edges_to_vtk(
                    const std::string & aFilePath,
                    uint aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            /**
             * Dumps the lagrange mesh into a VTK file
             *
             * @param[ in ] aFilePath  path of VTK file
             */
            void save_mesh_to_vtk(
                    const std::string & aFilePath,
                    uint aLagrangeMeshIndex );

            // -----------------------------------------------------------------------------

            void create_input_and_output_meshes();

            // -----------------------------------------------------------------------------

            void perform_initial_refinement();

            // -----------------------------------------------------------------------------

            bool get_restarted_from_file()
            {
                return mRestartedFromFile;
            };

        private:

            void user_defined_flagging(
                    Cell< hmr::Element * >   & aCells,
                    Cell< hmr::Element * >   & aCandidates,
                    const  Matrix< DDRMat >  & aVertexValues,
                    uint                       aFunctionIndex);

            void user_defined_flagging(
                    Cell< hmr::Element * >   & aCells,
                    Cell< hmr::Element * >   & aCandidates,
                    const  Matrix< DDRMat >  & aVertexValues,
                    Refinement_Function        aRefinementFunction);

    }; /* HMR */

} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_HPP_ */

