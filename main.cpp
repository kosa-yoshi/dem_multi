
#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "config.h"
#include "wall_boundary.h"
#include "dem.h"

int main(int argc, char **argv)
{
    Log::open("./log_simulation.txt", "dev_dem-calculation");

    // initialization
    ConfigParameter *config;
    config = new ConfigParameter();
    config->load_config_file();

    // developing
    DEM *dem;
    dem = new DEM();
    dem->initialize();
    dem->set_step(config->step);
    dem->open_log_file(config->calc_energy, config->colli_frequency, config->mes_contact_num, config->k_energy,
        config->calc_VanDerWaals_pp, config->calc_VanDerWaals_pw);

    Log::write("*** Simulation start ***\n");
    Log::write("calculation % .3e [s]\n", (double)config->total_step * dem->dt);
    clock_t sim_start  = clock();
    clock_t loop_start = clock();
    clock_t loop_end;

    while (true) {
        if (config->step > config->total_step) {
            break;
        }

        if (config->step >= config->output_step) {
            if (config->calc_dem)   {
                Log::write("write_ptcl_vtk\n");
                dem->write_dem_output(config->output_index);
            }

            if (config->export_vtk) {
                Log::write("write_wall_obj_vtk\n");
                dem->wall->write_vtk_file(config->output_index);
            }

            config->output_step  += config->output_frequency;
            config->output_index += 1;

            loop_end   = clock();
            Log::write("[output] calculation : step = % 8d/% 8d, time = % .2e [sec]\n",
                    config->step, config->total_step, (double)(loop_end - loop_start) / CLOCKS_PER_SEC);
            loop_start = clock();
        }


        dem->calculation();

        config->step += 1;
    }

    clock_t sim_end   = clock();
    Log::write("*** Simulation finished successfully !! ***\n");
    Log::write("Total calculation time : % .2e [sec]\n",
            (double)(sim_end - sim_start) / CLOCKS_PER_SEC);


    dem->finalization();
    delete dem;

    // finalization
    delete config;

    Log::close();
    return 0;
}


