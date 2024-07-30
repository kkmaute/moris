/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Reinitialize_Performer.hpp
 *
 */

#ifndef SRC_cl_WRK_Reinitialize_Performer
#define SRC_cl_WRK_Reinitialize_Performer

#include <memory>
#include "cl_Parameter_List.hpp"
#include "cl_Matrix.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace hmr
    {
        class HMR;
        class Mesh;
    }    // namespace hmr

    namespace mtk
    {
        class Field;
        class Mesh_Manager;
        class Mesh;
    }    // namespace mtk

    namespace MSI
    {
        enum class Dof_Type;
    }

    namespace gen
    {
        class Geometry_Engine;
    }

    namespace mdl
    {
        class Model;
    }

    namespace wrk
    {
        class Reinitialize_Performer
        {
          private:
            std::shared_ptr< Library_IO > mLibrary = nullptr;

            // Solution Dof types
            Vector< moris::MSI::Dof_Type > mDofTypes;

            // adv filed name
            std::string mADVFiledName;

            // discretization mesh index for adof
            moris_index mAdofMeshIndex;

            // coefficents that will be populated with target and source field coeff
            moris::Matrix< moris::DDRMat > mCoefficients;

            // reinitialization frequency
            sint mReinitializationFrequency;

            // mtk fields that will be used to overwrite gen fields
            Vector< std::shared_ptr< mtk::Field > > mMTKFields;

            // mesh outputting params
            std::string mOutputMeshFile;
            moris::real mTimeOffset;

          public:
            //------------------------------------------------------------------------------

            /**
             * @brief Construct a new Reinitialize_Performer object
             *
             * @param aParameterlist
             * @param aLibrary
             */

            Reinitialize_Performer( std::shared_ptr< Library_IO > aLibrary );

            //------------------------------------------------------------------------------

            /**
             * @brief Destroy the Reinitialize_Performer object
             *
             */

            ~Reinitialize_Performer(){};

            //------------------------------------------------------------------------------

            /**
             * @brief
             *
             * @param aHMRPerformers
             * @param aGENPerformer
             * @param aMTKPerformer
             * @param mMDLPerformer
             */

            void perform(
                    Vector< std::shared_ptr< hmr::HMR > >&            aHMRPerformers,
                    Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
                    Vector< std::shared_ptr< mtk::Mesh_Manager > >&   aMTKPerformer,
                    Vector< std::shared_ptr< mdl::Model > >           mMDLPerformer );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the reinitialization frequency object
             *
             * @return moris::uint
             */

            moris::sint
            get_reinitialization_frequency() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the coefficients object
             *
             * @return Matrix <DDRMat> const&
             */

            Matrix< DDRMat > const &
            get_coefficients() const;

            //------------------------------------------------------------------------------

            /**
             * @brief clip values of advs
             *
             * @param aGENPerformer
             */

            void
            impose_upper_lower_bound( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer, mtk::Field* aField );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the mtk fields object get the replaced fields to use in GEN
             *
             * @return Vector< std::shared_ptr< mtk::Field > >
             */
            Vector< std::shared_ptr< mtk::Field > >
            get_mtk_fields() const;

            //------------------------------------------------------------------------------

            /**
             * @brief  output the mapped and original field on the same mesh
             *
             * @param aTarget
             * @param aSource
             */
            void
            output_fields( mtk::Field* aTarget, mtk::Field* aSource, std::string aExoFileName ) const;
        };
    }    // namespace wrk
}    // namespace moris

#endif /* cl_WRK_Reinitailize_Performer.hpp */
