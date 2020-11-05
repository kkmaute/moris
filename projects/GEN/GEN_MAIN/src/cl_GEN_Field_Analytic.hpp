#ifndef MORIS_CL_GEN_FIELD_ANALYTIC_HPP
#define MORIS_CL_GEN_FIELD_ANALYTIC_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Analytic : virtual public Field
        {
            protected:
                bool mInterpolateChildNodes = false;
            public:

                /**
                 * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
                 */
                Field_Analytic();

                /**
                 * Given a node index or coordinate, returns the field value.
                 *
                 * @param aNodeIndex Node index
                 * @param aCoordinates Vector of coordinate values
                 * @return Field value
                 */
                real get_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

                /**
                 * Given a node coordinate, returns the field value
                 *
                 * @param aCoordinates vector of coordinate values
                 * @return Field value
                 */
                virtual real get_field_value_geometry(uint aNodeIndex,const Matrix<DDRMat>& aCoordinates) = 0;

                /**
                 * Given a node index or coordinate, returns a matrix of all sensitivities.
                 *
                 * @param aNodeIndex Node index
                 * @param aCoordinates Vector of coordinate values
                 * @return Matrix of sensitivities
                 */
                Matrix<DDRMat> get_field_sensitivities(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

                /**
                 * Given a node coordinate, returns a matrix of all sensitivities
                 *
                 * @param aCoordinates Vector of coordinate values
                 * @return Matrix of sensitivities
                 */
                virtual Matrix<DDRMat> get_field_sensitivities(const Matrix<DDRMat>& aCoordinates) = 0;

                void add_child_node(uint aNodeIndex, std::shared_ptr<Child_Node> aChildNode);

                void reset_child_nodes();

        };
    }
}

#endif //MORIS_CL_GEN_FIELD_ANALYTIC_HPP
