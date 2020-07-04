#ifndef DEM_CFD_PROJECT_CONFIG_H
#define DEM_CFD_PROJECT_CONFIG_H

#include "utility.h"
// #include "wall_boundary.h"

class ConfigParameter
{
    public:
        int step;
        int total_step;

        int output_frequency;
        int output_index;
        int output_step;

        int colli_frequency;

        bool calc_dem;
        bool calc_cfd;

        bool export_vtk;
        bool calc_energy;
        bool mes_contact_num;
        bool k_energy;
        bool calc_VanDerWaals_pp;
        bool calc_VanDerWaals_pw;

        ConfigParameter();
        ~ConfigParameter() {}

        void load_config_file();
};

#endif
