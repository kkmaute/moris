/*
 * cl_MTK_Vertex_Interpolation_DataBase.hpp
 *
 *  Created on: Dec  11, 2021
 *      Author: momo
 */
#ifndef SRC_cl_MTK_Vertex_Interpolation_DataBase_HPP_
#define SRC_cl_MTK_Vertex_Interpolation_DataBase_HPP_

#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;

        class Vertex_Interpolation_DataBase : public mtk::Vertex_Interpolation
        {
          private:
            moris_index mVertexIndex;
            moris_index mVertexOrder;
            mtk::Mesh*  mMesh = nullptr;

          public:
            // default constrctor
            Vertex_Interpolation_DataBase() = default;

            //------------------------------------------------------------------------------

            Vertex_Interpolation_DataBase( moris_index aVertexIndex,
                moris_index                            aVertexOrder,
                mtk::Mesh*                             aMesh );

            //------------------------------------------------------------------------------

            ~Vertex_Interpolation_DataBase() = default;

            //------------------------------------------------------------------------------
            /**
             * returns the IDs of the interpolation coefficients
             */
            virtual 
            Matrix< IdMat >
            get_ids() const override;

            //------------------------------------------------------------------------------

            /**
             * returns the indices of the interpolation coefficients
             */
            virtual
            Matrix< IndexMat >
            get_indices() const override;

            //------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this vertex
             */
            virtual 
            Matrix< IdMat >
            get_owners() const override;

            //------------------------------------------------------------------------------
            /**
             * set the interpolation weights
             */
            virtual 
            void
            set_weights( const moris::Matrix< DDRMat >& aWeights ) override;

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation weights
             */
            virtual 
            const Matrix< DDRMat >*
            get_weights() const override;


            //------------------------------------------------------------------------------
            /**
             * set the coefficient objects
             */
            virtual void
            set_coefficients( moris::Cell< Vertex* >& aCoefficients ) override;

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects
             */
            virtual 
            moris::Cell< Vertex* >&
            get_coefficients() override;

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects (const version)
             */
            virtual 
            const moris::Cell< Vertex* >&
            get_coefficients() const override;

            //------------------------------------------------------------------------------

            /**
             * returns the number of coefficients attributed to this basis
             */
            virtual 
            uint
            get_number_of_coefficients() const override;

            //------------------------------------------------------------------------------

            /**
             * @brief populate the data that is being passed as a refernce
             * 
             */

            void
            set_outward_data();

            //------------------------------------------------------------------------------

            /**
             * @brief memory usage of the vertex interpolations
             * 
             * @return size_t 
             */

            size_t
            capacity();
        };
    }// namespace mtk
}// namespace moris


#endif /* cl_MTK_Vertex_Interpolation_DataBase.hpp */