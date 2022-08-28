/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Stab_Penalty_Contact.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"     //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_Stab_Penalty_Contact : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                real mMasterVolume     = 0.5; // volume on master
                real mSlaveVolume      = 0.5; // volume on slave
                real mInterfaceSurface = 1.0; // surface on master/slave interface

                enum class SP_Property_Type
                {
                        MATERIAL,
                        MAX_ENUM
                };

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Stab_Penalty_Contact();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Stab_Penalty_Contact(){};

                //------------------------------------------------------------------------------
                /**
                 * reset the cluster measures required for this SP
                 */
                void reset_cluster_measures()
                {
                    // evaluate cluster measures from the cluster
                    mMasterVolume     = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::INTERP, mtk::Master_Slave::MASTER );
                    mSlaveVolume      = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::INTERP, mtk::Master_Slave::SLAVE );
                    mInterfaceSurface = mCluster->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Master_Slave::MASTER );
                    //                std::cout<<"mInterfaceSurface "<<mInterfaceSurface<<std::endl;
                    //                std::cout<<"mInterfaceSurface in FEM "<<mCluster->compute_volume()<<std::endl;
                }

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER )
                {
                    Stabilization_Parameter::set_dof_type_list( aDofTypes, aIsMaster );
                }

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                        moris::Cell< std::string >             & aDvStrings,
                        mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsMaster );
                }

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter value
                 */
                void eval_SP();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a master dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdMasterDOF ( 1 x numDerDof )
                 */
                void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdMasterDOF: not implemented." );
                }
                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a slave dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdSlaveDOF ( 1 x numDerDof )
                 */
                void eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdSlaveDOF: not implemented." );
                }
                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dPPdMasterDV ( 1 x numDerDv )
                 */
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdMasterDV: not implemented." );
                }

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a slave dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dSPdSlaveDV ( 1 x numDerDv )
                 */
                void eval_dSPdSlaveDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdSlaveDV: not implemented." );
                }
                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_ */

