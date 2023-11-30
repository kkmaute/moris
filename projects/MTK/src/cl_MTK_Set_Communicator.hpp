/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Set_Communicator.hpp
 *
 */
#pragma once

#include "containers.hpp"

namespace moris
{
    namespace mtk
    {
        // forward declare the mtk::Set
        class Set;
        class Block_Set;
        class Side_Set;
        class Double_Side_Set;

        //------------------------------------------------------------------------------

        /**
         * @brief This class handles communication of basic set information (e.g. element type, interpolation order, etc.)
         * for a list of sets. It has access to protected information of the mtk::Set's. This class is needed to perform
         * communication of set information in bulk rather than a set-by-set basis to improve parallel performance.
         */
        class Set_Communicator
        {
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            explicit Set_Communicator( moris::Cell< mtk::Set* >& aSetsToCommunicate );

            explicit Set_Communicator( moris::Cell< mtk::Block_Set* >& aBlocksToCommunicate )
                    : Set_Communicator( (moris::Cell< mtk::Set* >&)aBlocksToCommunicate ){};

            explicit Set_Communicator( moris::Cell< mtk::Side_Set* >& aSideSetsToCommunicate )
                    : Set_Communicator( (moris::Cell< mtk::Set* >&)aSideSetsToCommunicate ){};

            explicit Set_Communicator( moris::Cell< mtk::Double_Side_Set* >& aDoubleSideSetsToCommunicate )
                    : Set_Communicator( (moris::Cell< mtk::Set* >&)aDoubleSideSetsToCommunicate ){};

            //------------------------------------------------------------------------------

            ~Set_Communicator() {}

            //------------------------------------------------------------------------------

            bool check_sets( moris::Cell< mtk::Set* >& aSetsToCommunicate );

            //------------------------------------------------------------------------------
            void communicate( moris::Cell< mtk::Set* >& aSetsToCommunicate );
        };    // end class: mtk::Set_Communicator

        //------------------------------------------------------------------------------

    }    // namespace mtk

}    // namespace moris