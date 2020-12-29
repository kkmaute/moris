#ifndef MORIS_CL_GEN_VOXEL_HPP
#define MORIS_CL_GEN_VOXEL_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        class Voxel_Input : public Geometry, public Field_Analytic
        {

        private:

                moris::Matrix< DDRMat > mDomainDimensions;
                moris::Matrix< DDRMat > mDomainOffset;

                moris::Matrix< DDUMat > mVoxelField;
                moris::uint             mVoxelsInX;
                moris::uint             mVoxelsInY;
                moris::uint             mVoxelsInZ;

                moris::uint             mNumGrainInd;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Voxel_Input(
                    Matrix<DDRMat>&  aADVs,
                    Matrix<DDUMat>   aGeometryVariableIndices,
                    Matrix<DDUMat>   aADVIndices,
                    Matrix<DDRMat>   aConstants,
                    std::string      aVoxelFieldName,
                    Matrix<DDRMat>   aDomainDimensions,
                    Matrix<DDRMat>   aDomainOffset,
                    Field_Parameters aParameters = {});

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aFieldVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Voxel_Input(
                    sol::Dist_Vector* aOwnedADVs,
                    Matrix<DDUMat>    aGeometryVariableIndices,
                    Matrix<DDUMat>    aADVIndices,
                    Matrix<DDRMat>    aConstants,
                    std::string       aVoxelFieldName,
                    Matrix<DDRMat>    aDomainDimensions,
                    Matrix<DDRMat>    aDomainOffset,
                    Field_Parameters  aParameters = {});

            /**
             * Constructor with only constant parameters
             *
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Voxel_Input(
                    Matrix<DDRMat>   aConstants,
                    std::string      aVoxelFieldName,
                    Matrix<DDRMat>   aDomainDimensions,
                    Matrix<DDRMat>   aDomainOffset,
                    Field_Parameters aParameters = {});

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Distance to this geometry
             */
            real get_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

            moris::uint get_num_voxel_Ids(){ return mNumGrainInd; };

        private:

            void read_voxel_data( std::string aVoxelFieldName);


        };
    }
}

#endif //MORIS_CL_GEN_VOXEL_HPP
