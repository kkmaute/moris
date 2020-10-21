#ifndef MORIS_CL_GEN_Field_HPP
#define MORIS_CL_GEN_Field_HPP

#include "cl_Matrix.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace ge
    {
        class Field
        {
        protected:
            Cell<real*> mFieldVariables;

        private:
            sol::Dist_Vector* mSharedADVs = nullptr;
            Matrix<DDRMat> mConstantParameters;
            Matrix<DDSMat> mDeterminingADVIds;
            bool mDependsOnADVs;
            uint mNumADVs;
            sint mNumRefinements;
            sint mRefinementFunctionIndex;
            sint mBSplineMeshIndex;
            real mBSplineLowerBound;
            real mBSplineUpperBound;

        protected:

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full ADVs
             * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the field variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for B-spline discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(Matrix<DDRMat>& aADVs,
                  Matrix<DDUMat>  aFieldVariableIndices,
                  Matrix<DDUMat>  aADVIndices,
                  Matrix<DDRMat>  aConstantParameters,
                  sint            aNumRefinements,
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
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for B-spline discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(sol::Dist_Vector* aOwnedADVs,
                  Matrix<DDUMat>    aFieldVariableIndices,
                  Matrix<DDUMat>    aADVIndices,
                  Matrix<DDRMat>    aConstantParameters,
                  sint              aNumRefinements,
                  sint              aRefinementFunctionIndex,
                  sint              aBSplineMeshIndex,
                  real              aBSplineLowerBound,
                  real              aBSplineUpperBound);

            /**
             * Constructor for setting all field variables as consecutive ADVs.
             *
             * @param aSharedADVIds Shared ADV IDs needed for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for B-spline discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(const Matrix<DDSMat>& aSharedADVIds,
                  sint                  aNumRefinements,
                  sint                  aRefinementFunctionIndex,
                  sint                  aBSplineMeshIndex,
                  real                  aBSplineLowerBound,
                  real                  aBSplineUpperBound);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this field
             * @param aNumRefinements The number of refinement steps to use for this field
             * @param aRefinementFunctionIndex The index of a user-defined refinement function (-1 = default refinement)
             * @param aBSplineMeshIndex The index of a B-spline mesh for B-spline discretization (-1 = no B-splines)
             * @param aBSplineLowerBound The lower bound for the B-spline coefficients describing this field
             * @param aBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            Field(Matrix<DDRMat> aConstantParameters,
                  sint           aNumRefinements,
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
             * Resets all child nodes, called when a new XTK mesh is being created.
             */
            virtual void reset_child_nodes();

            /**
             * If this field depends on ADVs
             *
             * @return if this field has ADV indices
             */
            bool depends_on_advs();

            /**
             * Gets the IDs of ADVs which this field depends on for evaluations.
             *
             * @return Field basis ADV IDs
             */
            Matrix<DDSMat> get_determining_adv_ids();

            /**
             * This function will return true when called less than the number of refinements set for this field,
             * and false otherwise. This is to determine for a given refinement call if this field needs refinement.
             *
             * @return if to perform an additional refinement with this field
             */
            sint get_num_refinements();

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

            /**
             * Function for determining if this field is to be used for seeding a B-spline field.
             *
             * @return Logic for B-spline creation
             */
            virtual bool conversion_to_bsplines();

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
