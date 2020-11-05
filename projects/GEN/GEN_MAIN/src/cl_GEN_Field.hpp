#ifndef MORIS_CL_GEN_Field_HPP
#define MORIS_CL_GEN_Field_HPP

#include "cl_Matrix.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh;
    }
    namespace ge
    {
        class Field
        {
        protected:
            Cell<real*> mFieldVariables;

            Cell<std::shared_ptr<Child_Node>> mChildNodes;

            uint mNumOriginalNodes;
            mtk::Interpolation_Mesh* mMesh = nullptr;

        private:

            sol::Dist_Vector* mSharedADVs = nullptr;
            Matrix<DDRMat>    mConstantParameters;
            Matrix<DDSMat>    mDeterminingADVIds;
            std::string       mName;
            bool              mDependsOnADVs;
            Matrix<DDSMat>    mNumRefinements;
            Matrix<DDSMat>    mRefinementMeshIndices;
            sint              mRefinementFunctionIndex;
            sint              mBSplineMeshIndex;
            real              mBSplineLowerBound;
            real              mBSplineUpperBound;


        protected:

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full ADVs
             * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the field variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(Matrix<DDRMat>& aADVs,
                  Matrix<DDUMat>  aFieldVariableIndices,
                  Matrix<DDUMat>  aADVIndices,
                  Matrix<DDRMat>  aConstantParameters,
                  std::string     aName,
                  Matrix<DDSMat>  aNumRefinements,
                  Matrix<DDSMat>  aRefinementMeshIndices,
                  sint            aRefinementFunctionIndex,
                  sint            aBSplineMeshIndex,
                  real            aBSplineLowerBound,
                  real            aBSplineUpperBound);

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the field variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(sol::Dist_Vector* aOwnedADVs,
                  Matrix<DDUMat>    aFieldVariableIndices,
                  Matrix<DDUMat>    aADVIndices,
                  Matrix<DDRMat>    aConstantParameters,
                  std::string       aName,
                  Matrix<DDSMat>  aNumRefinements,
                  Matrix<DDSMat>  aRefinementMeshIndices,
                  sint              aRefinementFunctionIndex,
                  sint              aBSplineMeshIndex,
                  real              aBSplineLowerBound,
                  real              aBSplineUpperBound);

            /**
             * Constructor for setting all field variables as consecutive ADVs.
             *
             * @param aSharedADVIds Shared ADV IDs needed for this field
             * @param aName Name of this field for identification
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */

            Field(const Matrix<DDSMat>& aSharedADVIds,
                  std::string           aName,
                  Matrix<DDSMat>  aNumRefinements,
                  Matrix<DDSMat>  aRefinementMeshIndices,
                  sint                  aRefinementFunctionIndex,
                  sint                  aBSplineMeshIndex,
                  real                  aBSplineLowerBound,
                  real                  aBSplineUpperBound);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this field
             * @param aName Name of this field for identification
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(Matrix<DDRMat> aConstantParameters,
                  std::string    aName,
                  Matrix<DDSMat>  aNumRefinements,
                  Matrix<DDSMat>  aRefinementMeshIndices,
                  sint           aRefinementFunctionIndex,
                  sint           aBSplineMeshIndex,
                  real           aBSplineLowerBound,
                  real           aBSplineUpperBound);

            /**
             * Default constructor
             *
             * @error This only exists so the Blanca compiler does not complain, <b>DO NOT USE</b>.
             */
            Field();

        public:

            /**
             * Destructor
             */
            ~Field();

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            virtual real get_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates) = 0;

            /**
             * Given a node index or coordinate, returns a matrix of all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            virtual Matrix<DDRMat> get_field_sensitivities(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates) = 0;

            /**
             * Imports the local ADVs required from the full owned ADV distributed vector.
             *
             * @param aOwnedADVs Full owned distributed ADV vector
             */
            virtual void import_advs(sol::Dist_Vector* aOwnedADVs);

            /**
             * Add a new child node for evaluation, implemented for discrete fields
             *
             * @param aNodeIndex Index of the child node
             * @param aChildNode Contains information about how the child node was created
             */
            virtual void add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode);

            /**
             * Resets all nodal information, called when a new XTK mesh is being created.
             */
            virtual void reset_nodal_information();

            /**
             * Gets if this field is to be turned into a stored geometry/property, in order to store field values.
             *
             * @return Logic for storing field values
             */
            bool store_field_values();

            /**
             * Gets if this field is to be used for seeding a B-spline field.
             *
             * @return Logic for B-spline creation
             */
            virtual bool conversion_to_bsplines();

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Node coordinates
             * @return Determining ADV IDs at this node
             */
            virtual Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * If this field depends on ADVs
             *
             * @return if this field has ADV indices
             */
            bool depends_on_advs();

            /**
             * Gets the name of this field.
             *
             * @return Field name
             */
            std::string get_name();

            /**
             * This function will return true when called less than the number of refinements set for this field,
             * and false otherwise. This is to determine for a given refinement call if this field needs refinement.
             *
             * @return if to perform an additional refinement with this field
             */
            const Matrix< DDSMat > & get_num_refinements();

            const Matrix< DDSMat > & get_refinement_mesh_indices();

            /**
             * Gets the index of a user-defined refinement function used within HMR.
             *
             * @return User-defined refinement function index
             */
            sint get_refinement_function_index();

            /**
             * Gets the index of a B-Spline mesh for creating a B-Spline discretization for this field
             *
             * @return B-Spline mesh index, or -1 if not creating level set
             */
            sint get_bspline_mesh_index();

            /**
             * Gets the lower bound for the B-spline field.
             *
             * @return Lower bound
             */
            real get_bspline_lower_bound();

            /**
             * Get the upper bound for the B-spline field.
             *
             * @return Upper bound
             */
            real get_bspline_upper_bound();

            void set_mesh( mtk::Interpolation_Mesh* aMesh);

        private:

            /**
             * Checks variable inputs and resizes the internal field variables based these inputs.
             *
             * @param aFieldVariableIndices Indices of Field variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the Field variables
             */
            void assign_adv_dependencies(
                    Matrix<DDUMat>& aFieldVariableIndices,
                    Matrix<DDUMat>& aADVIndices);
            
            /**
             * Fills the remaining field variables with constant parameters.
             */
            void fill_constant_parameters();

        };
    }
}


#endif //MORIS_CL_GEN_Field_HPP
