#include "config.h"

ConfigParameter::ConfigParameter()
{
    step       = 0;
    total_step = 0;

    output_frequency = 0;
    output_index     = 0;
    output_step      = 0;

    colli_frequency = 0;

    calc_dem = false;
    calc_cfd = false;

    export_vtk  = false;
    calc_energy = false;
    mes_contact_num = false;
    k_energy        = false;
    calc_VanDerWaals_pp = false;
    calc_VanDerWaals_pw = false;
}

void ConfigParameter::load_config_file()
{
    FILE *fp_cfg;
    fp_cfg = Util::file_check_open("./input/config.txt", "r");

    Log::write("[INPUT] load config file !! (./input/config.txt)\n");

    Util::read_line_data(fp_cfg, "# calculation step");
    Util::get_int_data(fp_cfg, "start_index",      step);
    Util::get_int_data(fp_cfg, "total_steps",      total_step);

    Util::read_line_data(fp_cfg, "# output file");
    Util::get_int_data(fp_cfg, "frequency", output_frequency);
    Util::get_int_data(fp_cfg, "index",     output_index);

    Util::read_line_data(fp_cfg, "# config parameter");
    Util::get_bool_data(fp_cfg, "calc_dem", calc_dem);
    Util::get_bool_data(fp_cfg, "calc_cfd", calc_cfd);

    Util::get_bool_data(fp_cfg, "export_vtk",  export_vtk);
    Util::get_bool_data(fp_cfg, "calc_energy", calc_energy);
    Util::get_int_data( fp_cfg, "frequency_for_energy", colli_frequency);
    
    Util::get_bool_data(fp_cfg, "mes_contact_num", mes_contact_num);
    Util::get_bool_data(fp_cfg, "kinetic_energy", k_energy);
    Util::get_bool_data(fp_cfg, "VanDerWaals_pp", calc_VanDerWaals_pp);
    Util::get_bool_data(fp_cfg, "VanDerWaals_pw", calc_VanDerWaals_pw);

    fclose(fp_cfg);

    output_step = step;
}

