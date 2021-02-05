/*
 * cl_MTK_Field.hpp
 *
 *  Created on: Jan 19, 2021
 *      Author: schmidt
 */

#ifndef CL_MTK_FIELD_HPP_
#define CL_MTK_FIELD_HPP_

#include <memory>

#include "cl_Matrix.hpp"
#include "typedefs.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;
        class Mesh_Manager;

        class Field
        {
        private:
            uint mNumberOfDimensions = 1; //! Number of dimensions. Right now only 1.
            sint mDiscretizationMeshIndex; //! Discretization mesh index
            std::string mName; //! Name
            Matrix<DDRMat> mNodalValues;

        protected:
            /**
             * Constructor
             *
             * @param aDiscretizationMeshIndex Discretization mesh index (-1 if no discretization)
             * @param aName Name (optional)
             */
            Field(sint        aDiscretizationMeshIndex = -1,
                  std::string aName = "");

        public:

            /**
             * Destructor
             */
            virtual ~Field();

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
             * Gets all of the nodal values of this field.
             *
             * @param aMesh Mesh
             * @return Nodal values
             */
            const Matrix<DDRMat>& get_nodal_values(mtk::Mesh* aMesh);

            /**
             * Gets the number of dimensions of this field.
             *
             * @return Number of dimensions
             */
            uint get_number_of_dimensions();

            /**
             * Gets the name of this field.
             *
             * @return Name
             */
            std::string get_name();

            /**
             * Gets the index of a mesh where this field's discretization is defined.
             *
             * @return Mesh index (undefined behavior if discretization is not being used)
             */
            uint get_discretization_mesh_index();

        };
    }
}

#endif /* CL_MTK_FIELD_HPP_ */
