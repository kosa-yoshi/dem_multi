#include "dem.h"

DEM::DEM()
{
    wall = nullptr;
}

DEM::~DEM()
{
    delete wall;
}

void DEM::initialize()
{
    wall = new WallManager();
    wall->initialize();

    Log::write("[input] read dem config file (./input/config_dem.txt)\n");
    FILE *fp_dem;
    fp_dem = Util::file_check_open("./input/config_dem.txt", "r");

    Util::read_line_data( fp_dem, "# common parameter");
    Util::get_double_data(fp_dem, "time_step",   dt);
    Util::get_double_data(fp_dem, "spring",      spring);
    Util::get_double_data(fp_dem, "restitution", restitution);
    Util::get_double_data(fp_dem, "p_friction",  friction);
    Util::get_double_data(fp_dem, "w_friction",  friction_wall);
    Util::get_vec_data(   fp_dem, "gravity",     gravity);

    Util::get_char_data(  fp_dem, "div_axis",    axis);
    Util::get_int_data(   fp_dem, "num_species", num_species);
 
    num_particle_array = new int[num_species];
    diameter = new double[num_species];
    density  = new double[num_species];

    for (int i = 0; i < num_species; ++i)
    {
        Util::read_line_data( fp_dem, "# particle parameter");
        Util::get_int_data(   fp_dem, "particle_number",   num_particle_array[i]);
        Util::get_double_data(fp_dem, "particle_diameter", diameter[i]);
        Util::get_double_data(fp_dem, "particle_density",  density[i]);

        num_particle += num_particle_array[i];

        if (i != 0 && diameter[i - 1] > diameter[i])
        {
            for (int j = 0; j < i; ++j)
            {
                if (diameter[i] >= diameter[j])
                {
                    continue;
                }

                for (int k = j; k < i; ++k)
                {
                    double change = diameter[k];
                    diameter[k] = diameter[i];
                    diameter[i] = change;

                    int change_num = num_particle_array[k];
                    num_particle_array[k] = num_particle_array[i];
                    num_particle_array[i] = change_num;

                    change = density[k];
                    density[k] = density[i];
                    density[i] = change;
                }

                break;
            }
        }
    }

    Util::get_double_data(fp_dem, "Hamaker_pp", hamaker_pp);
    Util::get_double_data(fp_dem, "Hamaker_pw", hamaker_pw);
    Util::get_double_data(fp_dem, "cutoff_distance_min", cutoff_dist_min);
    Util::get_double_data(fp_dem, "cutoff_distance_max", cutoff_dist_max);

    Util::read_line_data(fp_dem, "# periodic boundary parameter");
    Util::get_double_data(fp_dem, "periodic_max", periodic_max);
    Util::get_double_data(fp_dem, "periodic_min", periodic_min);

    periodic_dist = fabs(periodic_max - periodic_min);
    Log::write("periodic_dist = % .4e\n", periodic_dist);
    
    // memory allocation
    pos    = new Vector3d[num_particle];
    vel    = new Vector3d[num_particle];
    force  = new Vector3d[num_particle];
    omega  = new Vector3d[num_particle];
    torque = new Vector3d[num_particle];

    species_cr = new int[num_particle];

    radius  = new double[num_species];
    volume  = new double[num_species];
    mass    = new double[num_species];
    moment  = new double[num_species];

    wall_contact = new int[num_particle];
    wall_strain  = new Vector3d[num_particle];

    for (int i = 0; i < num_particle; ++i) {
        wall_contact[i] = -1;  // -1 = no wall contact
    }

    // read particle data
    Util::read_line_data(fp_dem, "# initial particle data");
    int pid = 0;
    while (true) {
        char buf[256];
        if (!fgets(buf, sizeof(buf), fp_dem)) {
            Log::write("*** finish reading initial particle data\n");
            break;
        }

        Vector3d init_pos;
        Vector3d init_vel;
        Vector3d init_rot;

        int num_scaned = sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &init_pos.x, &init_pos.y, &init_pos.z,
                &init_vel.x, &init_vel.y, &init_vel.z,
                &init_rot.x, &init_rot.y, &init_rot.z);

        if (pid < num_particle) {
            if (num_scaned != 9) {
                if (num_scaned == 3) {
                    init_vel.initialize();
                    init_rot.initialize();
                } else if (num_scaned == 6) {
                    init_rot.initialize();
                } else {
                    Log::write("[ERROR] unexpected buf data line# number of scaned data : [%d]\n", num_scaned);
                    exit(1);
                }
            }

            pos[pid]   = init_pos;
            vel[pid]   = init_vel;
            omega[pid] = init_rot;

            ++pid;
        } else {
            Log::write("[Warning] too much initial data !!\ncannot read all particle data\n");
            break;
        }
    }

    fclose(fp_dem);

    // generate particle data
    Vector3d min_pos, max_pos;
    wall->get_min_point_position(min_pos);
    wall->get_max_point_position(max_pos);

    int num_generate_ptcl_sum = num_particle - pid;

    if (num_generate_ptcl_sum > 0) {
        Log::write("=========================\n");
        Log::write("** particle generation **\n");

        for (int j = 0; j < num_species; ++j)
        {
            int num_generate_ptcl = num_generate_ptcl_sum - num_particle_array[j];

            if (num_generate_ptcl > 0)
            {
                num_generate_ptcl = num_particle_array[j];
            }

            else
            {
                num_generate_ptcl = num_generate_ptcl_sum;
            }

            int offset = 1;
            int nx = (int)floor((max_pos.x - min_pos.x) / diameter[j]) - offset * 2;
            int ny = (int)floor((max_pos.y - min_pos.y) / diameter[j]) - offset * 2;
            int nz = (int)floor((max_pos.z - min_pos.z) / diameter[j]) - offset * 2;

            int num_max = nx * ny * nz;
            if (num_generate_ptcl > num_max) {
                Log::write("too much number of generate particles\n");
                Log::write("please set the total particle number under [%d]\n", num_max);
                exit(1);
            }

            Log::write("particle generate species No.%d\n", j);
            Log::write("particle generate from [%d] == ", pid);

            for (int id = 0; id < num_generate_ptcl; ++id) {
                //if (id == num_max) {
                  //  break;
                //}

                int px = id % nx;
                int py = id / (nx * nz);
                int pz = (id / nx) % nz;

                pos[pid].x = (0.5 + (double)(px + offset)) * diameter[j] + min_pos.x;
                pos[pid].y = (0.5 + (double)(py + offset)) * diameter[j] + min_pos.y;
                pos[pid].z = (0.5 + (double)(pz + offset)) * diameter[j] + min_pos.z;

                vel[pid].x = ((rand() % 2) == 0) ? (0.0001 * (rand() % 25)) : (-0.0001 * (rand() % 25));
                vel[pid].y = ((rand() % 2) == 0) ? (0.0001 * (rand() % 25)) : (-0.0001 * (rand() % 25));
                vel[pid].z = ((rand() % 2) == 0) ? (0.0001 * (rand() % 25)) : (-0.0001 * (rand() % 25));

                ++pid;

            }

            Log::write(" to [%d]\n", pid);
            Log::write("=========================\n");

            num_generate_ptcl_sum -= num_generate_ptcl;
        }
    }

    wall->set_object_param(dt);

    // initialize contact list
    const int num_possible_contact = 12;

    ptcl_contact_num = new int[num_particle];
    Util::allocate_int_array_2d(ptcl_contact_list, num_particle, num_possible_contact);
    Util::allocate_int_array_2d(ptcl_contact_list_prev, num_particle, num_possible_contact);
    Util::allocate_vector_array_2d(ptcl_contact_strain, num_particle, num_possible_contact);
    Util::allocate_vector_array_2d(ptcl_contact_strain_prev, num_particle, num_possible_contact);

    for (int i = 0; i < num_particle; ++i) {
        ptcl_contact_num[i] = 0;
        for (int j = 0; j < num_possible_contact; ++j) {
            ptcl_contact_list[i][j]      = -1;
            ptcl_contact_list_prev[i][j] = -1;

            ptcl_contact_strain[i][j].initialize();
            ptcl_contact_strain_prev[i][j].initialize();
        }
    }

    // initialize near particle list
    celldata = new Cell[num_species];

    for (int i = 0; i < num_species; ++i)
    {
        celldata[i].cell_size = diameter[i];

        celldata[i].num_cell_x = (int)((max_pos.x - min_pos.x) / celldata[i].cell_size) + 2;
        celldata[i].num_cell_y = (int)((max_pos.y - min_pos.y) / celldata[i].cell_size) + 2;
        celldata[i].num_cell_z = (int)((max_pos.z - min_pos.z) / celldata[i].cell_size) + 2;

        celldata[i].num_cell = celldata[i].num_cell_x * celldata[i].num_cell_y * celldata[i].num_cell_z;

        celldata[i].cell_pos.x = -min_pos.x + celldata[i].cell_size;
        celldata[i].cell_pos.y = -min_pos.y + celldata[i].cell_size;
        celldata[i].cell_pos.z = -min_pos.z + celldata[i].cell_size;

        celldata[i].cell_number    = new int[num_particle];
        celldata[i].count_per_cell = new int[celldata[i].num_cell];
        celldata[i].sort_list      = new int[num_particle];
    }

    // calculate physical properties
    for (int i = 0; i < num_species; ++i)
    {
        radius[i] = 0.5 * diameter[i];
        volume[i] = (3.14 / 6.0) * diameter[i] * diameter[i] * diameter[i];
        mass[i]   = volume[i] * density[i];
        moment[i] = 0.4 * mass[i] * radius[i] * radius[i];
    }

    int calc_sp = num_species + 1;
    Util::allocate_double_array_2d(converted_mass, num_species, calc_sp);
    Util::allocate_double_array_2d(converted_diameter, num_species, calc_sp);
    Util::allocate_double_array_2d(dashpot, num_species, calc_sp);

    for (int i = 0; i < num_species; ++i)
    {
        for (int j = 0; j < num_species; ++j)
        {
            converted_diameter[i][j] = diameter[i] * diameter[j] / (diameter[i] + diameter[j]);
            converted_mass[i][j] = mass[i] * mass[j] / (mass[i] + mass[j]);
        }

        // for wall
        converted_diameter[i][num_species] = diameter[i];
        converted_mass[i][num_species] = mass[i];
    }

    for (int i = 0; i < num_species; ++i)
    {
        for (int j = 0; j < calc_sp; ++j)
        {
            dashpot[i][j] = -2.0 * log(restitution) * sqrt((spring * converted_mass[i][j]) / (log(restitution) * log(restitution) + 3.14 * 3.14));
        }
    }

    for (int i = 0; i < num_species; ++i)
    {
        double limit_dt = 0.1 * 3.14 * sqrt(mass[i] / spring);
        if (dt > limit_dt) {
            Log::write("*** Warning !! time step (dt) should be set smaller ***\n");
            Log::write(" - dt       = % .4e\n", dt);
            Log::write(" - limit_dt = % .4e\n", limit_dt);
            // exit(1);  // for debug only
        }
    }

    for (int i = 0; i < num_species; ++i)
    {
        Log::write("ptcl.radius : % .4e\n", radius[i]);
        Log::write("ptcl.volume : % .4e\n", volume[i]);
        Log::write("ptcl.mass   : % .4e\n", mass[i]);
        Log::write("ptcl.moment : % .4e\n", moment[i]);
        //Log::write("ptcl.dashpot: % .4e\n", dashpot[i]);
    }
    /*
    // define each species position 
    int axis_flag;
    std::string str_axis = std::string(axis);

    if (str_axis == "y") {
        Log::write("sort along Y axis\n");
        axis_flag = 1;
    }
    else if (str_axis == "z") {
        Log::write("sort along Z axis\n");
        axis_flag = 2;
    }
    else {
        Log::write("sort along X axis\n");
        axis_flag = 0;
    }
  
    if (axis_flag == 0)
    {
        for (int j = 0; j < num_particle; ++j)
        {
            int addup = 0;
            for (int i = 0; i < num_particle; ++i)
            {
                if (i == j)
                {
                    continue;
                }

                if (pos[j].x > pos[i].x)
                {
                    addup += 1;
                }
                
                else if (pos[j].x == pos[i].x)
                {
                    if(j > i)
                    {
                        addup += 1;
                    }
                }
            }

            int num_particle_addup = 0;
            for (int n = 0; n < num_species; ++n)
            {
                num_particle_addup += num_particle_array[n];
                if (addup < num_particle_addup)
                {
                    species_cr[j] = n;
                    break;
                }
            }
        }
    }

    if (axis_flag == 1)
    {
        for (int j = 0; j < num_particle; ++j)
        {
            int addup = 0;
            for (int i = 0; i < num_particle; ++i)
            {
                if (i == j)
                {
                    continue;
                }

                if (pos[j].y > pos[i].y)
                {
                    addup += 1;
                }

                else if (pos[j].y == pos[i].y)
                {
                    if (j > i)
                    {
                        addup += 1;
                    }
                }
            }

            int num_particle_addup = 0;
            for (int n = 0; n < num_species; ++n)
            {
                num_particle_addup += num_particle_array[n];
                if (addup < num_particle_addup)
                {
                    species_cr[j] = n;
                    break;
                }
            }
        }
    }

    if (axis_flag == 2)
    {
        for (int j = 0; j < num_particle; ++j)
        {
            int addup = 0;
            for (int i = 0; i < num_particle; ++i)
            {
                if (i == j)
                {
                    continue;
                }

                if (pos[j].z > pos[i].z)
                {
                    addup += 1;
                }

                else if (pos[j].z == pos[i].z)
                {
                    if (j > i)
                    {
                        addup += 1;
                    }
                }
            }

            int num_particle_addup = 0;
            for (int n = 0; n < num_species; ++n)
            {
                num_particle_addup += num_particle_array[n];
                if (addup < num_particle_addup)
                {
                    species_cr[j] = n;
                    break;
                }
            }
        }
    }
    */

    for (int i = 0; i < num_particle; ++i)
    {
        if (i < num_particle_array[0])
        {
            species_cr[i] = 0;
        }

        else
        {
            species_cr[i] = 1;
        }
    }

    FILE* species_info;
    species_info = Util::file_check_open("separate_file_criterion.txt", "w");
    fprintf(species_info, "species_num:::::: %8d\n", num_species);
    for (int i = 0; i < num_species; ++i)
    {
        fprintf(species_info, "#particle_information\n");
        fprintf(species_info, "particle_num::::: %8d\n", num_particle_array[i]);
        fprintf(species_info, "radius::::::::::: %.4e\n", radius[i]);
        fprintf(species_info, "volume::::::::::: %.4e\n", volume[i]);
        fprintf(species_info, "mass::::::::::::: %.4e\n", mass[i]);
    }

    fprintf(species_info, "#species_information\n");
    for (int i = 0; i < num_particle; ++i)
    {
        fprintf(species_info, "species_cr::::::: %8d\n", species_cr[i]);
    }

    fclose(species_info);
}

void DEM::finalization()
{
    close_log_file();

    delete [] pos;
    delete [] vel;
    delete [] force;
    delete [] omega;
    delete [] torque;

    delete [] ptcl_contact_num;
    Util::free_int_array_2d(ptcl_contact_list);
    Util::free_int_array_2d(ptcl_contact_list_prev);
    Util::free_vector_array_2d(ptcl_contact_strain);
    Util::free_vector_array_2d(ptcl_contact_strain_prev);

    delete [] wall_contact;
    delete [] wall_strain;

    delete [] species_cr;
    delete [] num_particle_array;

    delete [] diameter;
    delete [] radius;
    delete [] volume;
    delete [] density;
    delete [] mass;
    delete [] moment;

    Util::free_double_array_2d(dashpot);
    Util::free_double_array_2d(converted_mass);
    Util::free_double_array_2d(converted_diameter);

    delete [] celldata;
}

void DEM::write_dem_output(const int output_id)
{
    FILE* vtk1;
    char file_name[64];
    sprintf(file_name, "./output/dem_ptcl%04d.vtk", output_id);
    vtk1 = Util::file_check_open(file_name, "w");

    fprintf(vtk1, "# vtk DataFile Version 3.0\n");
    fprintf(vtk1, "dem_particle_data\n");
    fprintf(vtk1, "ASCII\n");
    fprintf(vtk1, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(vtk1, "POINTS %d double\n", num_particle);
    for (int i = 0; i < num_particle; ++i) {
        fprintf(vtk1, "% .4e % .4e % .4e\n", pos[i].x, pos[i].y, pos[i].z);
    }
    fprintf(vtk1, "POINT_DATA %d\n", num_particle);

    fprintf(vtk1, "VECTORS vel double\n");
    for (int i = 0; i < num_particle; ++i) {
        fprintf(vtk1, "% .4e % .4e % .4e\n", vel[i].x, vel[i].y, vel[i].z);
    }

    fprintf(vtk1, "VECTORS force double\n");
    for (int i = 0; i < num_particle; ++i) {
        fprintf(vtk1, "% .4e % .4e % .4e\n", force[i].x, force[i].y, force[i].z);
    }

    fclose(vtk1);
    
    // export vtk data
    for (int j = 0; j < num_species; ++j)
    {
        FILE* vtk;
        //char file_name[64];
        sprintf(file_name, "./output/dem_ptcl%02d_%04d.vtk", j, output_id);
        vtk = Util::file_check_open(file_name, "w");

        fprintf(vtk, "# vtk DataFile Version 3.0\n");
        fprintf(vtk, "dem_particle_data\n");
        fprintf(vtk, "ASCII\n");
        fprintf(vtk, "DATASET UNSTRUCTURED_GRID\n");

        fprintf(vtk, "POINTS %d double\n", num_particle_array[j]);
        for (int i = 0; i < num_particle; ++i) {
            if (j == species_cr[i])
            {
                fprintf(vtk, "% .4e % .4e % .4e\n", pos[i].x, pos[i].y, pos[i].z);
            }
        }
        fprintf(vtk, "POINT_DATA %d\n", num_particle_array[j]);

        fprintf(vtk, "VECTORS vel double\n");
        for (int i = 0; i < num_particle; ++i) {
            if (j == species_cr[i])
            {
                fprintf(vtk, "% .4e % .4e % .4e\n", vel[i].x, vel[i].y, vel[i].z);
            }
        }

        fprintf(vtk, "VECTORS force double\n");
        for (int i = 0; i < num_particle; ++i) {
            if (j == species_cr[i])
            {
                fprintf(vtk, "% .4e % .4e % .4e\n", force[i].x, force[i].y, force[i].z);
            }
        }

        fclose(vtk);
    }

    // write restart file
    FILE *out;
    out = Util::file_check_open("./output/_restart_dem_config.txt", "w");

    fprintf(out, "# common parameter\n");
    fprintf(out, "time_step:::::::::: % .4e\n", dt);
    fprintf(out, "spring::::::::::::: % .4e\n", spring);
    fprintf(out, "restitution:::::::: % .4e\n", restitution);
    fprintf(out, "p_friction::::::::: % .4e\n", friction);
    fprintf(out, "w_friction::::::::: % .4e\n", friction_wall);
    fprintf(out, "gravity_(x_y_z):::: % .4e % .4e % .4e\n", gravity.x, gravity.y, gravity.z);

    fprintf(out, "div_axis::::::::::: %s\n", axis);
    fprintf(out, "num_species:::::::: %d\n", num_species);

    for (int i = 0; i < num_species; ++i)
    {
        fprintf(out, "# particle parameter\n");
        fprintf(out, "particle_number:::: %8d\n", num_particle_array[i]);
        fprintf(out, "particle_diameter:: % .4e\n", diameter[i]);
        fprintf(out, "particle_density::: % .4e\n", density[i]);
    }

    fprintf(out, "Hamaker_pp::::::::: % .4e\n", hamaker_pp);
    fprintf(out, "Hamaker_pw::::::::: % .4e\n", hamaker_pw);
    fprintf(out, "cutoff_distance_min::: % .4e\n", cutoff_dist_min);
    fprintf(out, "cutoff_distance_max::: % .4e\n", cutoff_dist_max);

    fprintf(out, "# periodic boundary parameter\n");
    fprintf(out, "periodic_max::::::: % .4e\n", periodic_max);
    fprintf(out, "periodic_min::::::: % .4e\n", periodic_min);

    fprintf(out, "# initial particle data\n");
    for (int i = 0; i < num_particle; ++i) {
        fprintf(out, "% .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e\n",
                pos[i].x, pos[i].y, pos[i].z, vel[i].x, vel[i].y, vel[i].z, omega[i].x, omega[i].y, omega[i].z);
    }

    fclose(out);
}

void DEM::open_log_file(const bool b_energy, const int colli_frequency, const bool b_contact_num, const bool b_kenergy, const bool b_van, const bool b_van_pw)
{
    k_energy = b_kenergy;
    mes_contact_num = b_contact_num;
    calc_energy = b_energy;
    calc_VanDerWaals_pp = b_van;
    calc_VanDerWaals_pw = b_van_pw;
    colli_output_freq = colli_frequency;

    if (calc_energy) {
        Log::write("open log file for collision energy\n");
        fp_energy_ptcl = Util::file_check_open("./output/_log_collision_energy_ptcl.csv", "w");
        fp_energy_obj  = Util::file_check_open("./output/_log_collision_energy_obj.csv",  "w");

        fp_colli_ptcl = Util::file_check_open("./output/_log_processed_collision_energy_ptcl.csv", "w");
        fp_colli_wall = Util::file_check_open("./output/_log_processed_collision_energy_obj.csv", "w");

        fprintf(fp_energy_ptcl, "step,i,j,pos_x,pos_y,pos_z,col_en,col_et\n");
        fprintf(fp_energy_obj,  "step,i,j,pos_x,pos_y,pos_z,col_en,col_et\n");
        
        fprintf(fp_colli_ptcl, "step, ave_col_en, max_col_en, min_col_en, ave_col_et, max_col_et, min_col_et\n");
        fprintf(fp_colli_wall, "step, ave_col_en, max_col_en, min_col_en, ave_col_et, max_col_et, min_col_et\n");
        
    }

    if (mes_contact_num) {
        Log::write("open log file for contact number\n");
        fp_contact_num_ptcl = Util::file_check_open("./output/_log_contact_num_ptcl.csv", "w");

        fprintf(fp_contact_num_ptcl, "step,contact_num_sum\n");
    }

    if(k_energy)
    {
        Log::write("open log file for kinetic energy\n");
        fp_kinetic_energy = Util::file_check_open("./output/_log_kinetic_energy.csv", "w");
        fprintf(fp_kinetic_energy, "step, average kinetic energy, max, min\n");
    }

    fp_over = Util::file_check_open("./output/_overlap.csv", "w");
    if (calc_VanDerWaals_pp)
    {
        Log::write("open log file for VanDerWaals force between particles\n");
    //    fp_van_pp = Util::file_check_open("./output/_log_VanDerWaals_pp.csv", "w");
    //    fprintf(fp_van_pp, "step, i, j, VanDerWaals force [x, y, z]\n");
    }

    if (calc_VanDerWaals_pw)
    {
        Log::write("open log file for VanDerWaals force between particle and wall\n");
    //    fp_van_pw = Util::file_check_open("./output/_log_VanDerWaals_pw.csv", "w");
    //    fprintf(fp_van_pw, "step, i, VanDerWaals force [x, y, z]\n");
    }
}

void DEM::close_log_file()
{
    if (calc_energy) {
        fclose(fp_energy_ptcl);
        fclose(fp_energy_obj);
        fclose(fp_colli_ptcl);
        fclose(fp_colli_wall);
    }

    if (mes_contact_num) {
        fclose(fp_contact_num_ptcl);
    }

    if (k_energy)
    {
        fclose(fp_kinetic_energy);
    }

    fclose(fp_over);

    //if (calc_VanDerWaals_pp)
    //{
    //    fclose(fp_van_pp);
    //}

    //if (calc_VanDerWaals_pw)
    //{
    //    fclose(fp_van_pw);
    //}
}

void DEM::build_near_particle_list()
{
    for (int j = 0; j < num_species; ++j)
    {
        double per_size = 1.0 / celldata[j].cell_size;

        for (int i = 0; i < num_particle; ++i) {
            int x = (int)floor((pos[i].x + celldata[j].cell_pos.x) * per_size);
            int y = (int)floor((pos[i].y + celldata[j].cell_pos.y) * per_size);
            int z = (int)floor((pos[i].z + celldata[j].cell_pos.z) * per_size);

            if (x <= 0) {
                x = 1;
            }

            if (x >= (celldata[j].num_cell_x - 1)) {
                x = celldata[j].num_cell_x - 2;
            }

            if (y <= 0) {
                y = 1;
            }

            if (y >= (celldata[j].num_cell_y - 1)) {
                y = celldata[j].num_cell_y - 2;
            }

            if (z <= 0) {
                z = 1;
            }

            if (z >= (celldata[j].num_cell_z - 1)) {
                z = celldata[j].num_cell_z - 2;
            }

            celldata[j].cell_number[i] = x + celldata[j].num_cell_x * y + celldata[j].num_cell_x * celldata[j].num_cell_y * z;
        }

        for (int i = 0; i < celldata[j].num_cell; ++i) {
            celldata[j].count_per_cell[i] = 0;
        }
       
        for (int i = 0; i < num_particle; ++i) {
            celldata[j].count_per_cell[celldata[j].cell_number[i]] += 1;
        }

        for (int i = 1; i < celldata[j].num_cell; ++i) {
            celldata[j].count_per_cell[i] += celldata[j].count_per_cell[i - 1];
        }

        for (int i = num_particle - 1; i >= 0; i--) {
            celldata[j].sort_list[--celldata[j].count_per_cell[celldata[j].cell_number[i]]] = i;
        }
    }
}

void DEM::calculate_ptcl_collision(const int i, const int j, int& contact, const double dist, const Vector3d& norm)
{
    Vector3d strain;
    bool is_first_collision = true;
    for (int n = 0; n < ptcl_contact_num[i]; ++n) {
        if (j == ptcl_contact_list_prev[i][n]) {
            strain = ptcl_contact_strain_prev[i][n];
            is_first_collision = false;
            break;
        }
    }

    // calculate obj velocity
    Vector3d rel_vel = vel[i] - vel[j];

    // normal direction
    double   norm_dot = norm * rel_vel;
    Vector3d rel_vel_n = norm_dot * norm;

    Vector3d overlap_n = (dist - radius[species_cr[i]] - radius[species_cr[j]]) * norm;
    Vector3d force_n = -spring * overlap_n - dashpot[species_cr[i]][species_cr[j]] * rel_vel_n;

    // tangential direction
    Vector3d rel_vel_t = rel_vel - rel_vel_n + ((radius[species_cr[i]] * (omega[i] + omega[j])) % norm);
    Vector3d overlap_t = rel_vel_t * dt + strain;
    Vector3d force_t = -spring * overlap_t - dashpot[species_cr[i]][species_cr[j]] * rel_vel_t;

    // check friction slider
    double abs_fn = friction * force_n.norm();
    double abs_ft = force_t.norm();
    if (abs_ft > abs_fn) {
        // if friction slider works, force in a tangential direction is calculated by a friction slider.
        Vector3d norm_t = rel_vel_t;
        double scale = norm_t.norm();
        if (scale == 0.0) {
            scale = strain.norm();
            // if (scale == 0.0) ERROR(unexpected error: rel_vel_t == strain == 0.0)
        }
        norm_t.x /= scale;
        norm_t.y /= scale;
        norm_t.z /= scale;

        force_t = -abs_fn * norm_t; // check
    }
    else {
        // if no friction slider works, overlap in a tangential direction is stored as a strain.
        strain = overlap_t;
    }

    // calculate force and torque
    Vector3d total_force = force_n + force_t;
    force[i] += total_force;
    force[j] -= total_force;

    Vector3d total_torque = radius[species_cr[i]] * (norm % force_t);
    torque[i] += total_torque;
    torque[j] += total_torque;

    // contact list update
    ptcl_contact_list[i][contact] = j;
    ptcl_contact_strain[i][contact] = strain;
    ++contact;

    double overlap_checki = overlap_n.norm();

    if (overlap_checki > 0.05 * radius[species_cr[i]] || overlap_checki > 0.05 * radius[species_cr[j]])
    {
        fprintf(fp_over, "i %d, % d, %d, %.4e\n", step, i, j, overlap_checki);
    }

    // output collision energy
    if (is_first_collision && calc_energy) {
        if (step % colli_output_freq == 0)
        {
            Vector3d center = 0.5 * (pos[i] + pos[j]);
            double colllision_en = 0.5 * converted_mass[species_cr[i]][species_cr[j]] * (rel_vel_n.x * rel_vel_n.x + rel_vel_n.y * rel_vel_n.y + rel_vel_n.z * rel_vel_n.z);
            double colllision_et = 0.5 * converted_mass[species_cr[i]][species_cr[j]] * (rel_vel_t.x * rel_vel_t.x + rel_vel_t.y * rel_vel_t.y + rel_vel_t.z * rel_vel_t.z);
            fprintf(fp_energy_ptcl, "%8d,%8d,%8d,% .4e,% .4e,% .4e,% .4e,% .4e\n",
                step, i, j, center.x, center.y, center.z, colllision_en, colllision_et);
        }
    }

    //output conatct number
    if (mes_contact_num)
    {
        if (species_cr[i] == 0 && species_cr[j] == 0)
        {
            count_00 += 1;
        }
        else if (species_cr[i] == 0 && species_cr[j] == 1)
        {
            count_01 += 1;
        }
        else if(species_cr[i] == 1 && species_cr[j] == 1)
        {
            count_11 += 1;
        }
    }
}

void DEM::print_output()
{
    if (mes_contact_num) 
    {
        fprintf(fp_contact_num_ptcl, "%8d,%8d,%8d,%8d\n", step, count_00, count_01, count_11);
        count_00 = 0;
        count_01 = 0;
        count_11 = 0;
    }

    if (k_energy)
    {
        ave_k_energy = kinetic_energy / num_particle;
        fprintf(fp_kinetic_energy, "%8d,%.4e,%.4e,%.4e\n", step, ave_k_energy, max_k_energy, min_k_energy);

        ave_k_energy = 0.0;
        kinetic_energy = 0.0;
    }
}

void DEM::calculate_ptcl_kinetic_energy()
{
    for (int i = 0; i < num_particle; ++i) {

        kinetic_energy += 0.5 * mass[species_cr[i]] * vel[i] * vel[i];
        kinetic_energy_pre = 0.5 * mass[species_cr[i]] * vel[i] * vel[i];

        if (i == 0)
        {
            min_k_energy = kinetic_energy_pre;
            max_k_energy = kinetic_energy_pre;
        }
        else
        {
            if (max_k_energy < kinetic_energy_pre)
            {
                max_k_energy = kinetic_energy_pre;
            }
            else if (min_k_energy > kinetic_energy_pre)
            {
                min_k_energy = kinetic_energy_pre;
            }
        }
    }
}

void DEM::calculate_VanDerWaals_pp(const int i, const int j, const double distance_bet_surface,const Vector3d& norm)
{
    // position vector
    //Vector3d direction_pp = norm; 
    //double converted_diameter = diameter[species_cr[i]] * diameter[species_cr[j]] / ( diameter[species_cr[i]] + diameter[species_cr[j]]);
    double vanderwaals_pp = hamaker_pp * converted_diameter[species_cr[i]][species_cr[j]] / 12 / distance_bet_surface / distance_bet_surface;

    Vector3d vanderwaals_pp_n = vanderwaals_pp * norm;

    force[i] -= vanderwaals_pp_n;
    force[j] += vanderwaals_pp_n;
    //fprintf(fp_van_pp, "%8d,%8d,%8d,%.4e,%.4e,%.4e,%.4e\n", step, i, j, vanderwaals_pp_n.x, vanderwaals_pp_n.y, vanderwaals_pp_n.z, distance_bet_surface);
}

void DEM::calculate_VanDerWaals_pw(const int i, const double distance_bet_surface, const Vector3d& wall_norm)
{
    // position vector
    //Vector3d direction_pw = wall_norm / wall_norm.norm();
    double vanderwaals_pw = hamaker_pw * diameter[species_cr[i]] / 12 / distance_bet_surface / distance_bet_surface;

    Vector3d vanderwaals_pw_n = vanderwaals_pw * wall_norm;

    force[i] -= vanderwaals_pw_n;
    //fprintf(fp_van_pw, "%8d,%8d,%.4e,%.4e,%.4e,%.4e\n", step, i, vanderwaals_pw_n.x, vanderwaals_pw_n.y, vanderwaals_pw_n.z, distance_bet_surface);
}


void DEM::calculate_wall_collision(const int i, const int obj_idx, const double dist, const Vector3d& wall_norm)
{
    Vector3d strain;
    bool is_first_collision = true;
    if (obj_idx == wall_contact[i]) {
        // collision at the same obj in previous time step
        strain = wall_strain[i];
        is_first_collision = false;
    } else {
        // collision at the other obj in previous time step or at the first time
        strain.initialize();
    }

    // calculate obj velocity
    Vector3d wall_vel;      // This code only support rotation velocity.
    Vector3d wall_pos = pos[i] + dist * wall_norm;
    if (wall->obj[obj_idx].num_axis > 0) {
        wall->obj[obj_idx].get_wall_velocity(wall_pos, wall_vel);
    }

    Vector3d rel_vel = vel[i] - wall_vel;

    // normal direction
    double   norm_dot  = wall_norm * rel_vel;
    Vector3d rel_vel_n = norm_dot  * wall_norm;

    Vector3d overlap_n = (dist - radius[species_cr[i]]) * wall_norm;
    Vector3d force_n   = -spring * overlap_n - dashpot[species_cr[i]][num_species] * rel_vel_n;

    // tangential direction
    Vector3d rel_vel_t = rel_vel - rel_vel_n + ((radius[species_cr[i]] * omega[i]) % wall_norm);
    Vector3d overlap_t = rel_vel_t * dt + strain;
    Vector3d force_t   = -spring * overlap_t - dashpot[species_cr[i]][num_species] * rel_vel_t;

    // check friction slider
    double abs_fn = friction_wall * force_n.norm();
    double abs_ft = force_t.norm();
    if (abs_ft > abs_fn) {
        // if friction slider works, force in a tangential direction is calculated by a friction slider.
        Vector3d norm_t = rel_vel_t;
        double scale = norm_t.norm();
        if (scale == 0.0) {
            scale = strain.norm();
            // if (scale == 0.0) ERROR(unexpected error: rel_vel_t == strain == 0.0)
        }
        norm_t.x /= scale;
        norm_t.y /= scale;
        norm_t.z /= scale;

        force_t = -abs_fn * norm_t; // check
    } else {
        // if no friction slider works, overlap in a tangential direction is stored as a strain.
        strain = overlap_t;
    }

    // calculate force and torque
    force[i]  += force_n + force_t;
    torque[i] += radius[species_cr[i]] * (wall_norm % force_t);

    // contact list update
    wall_contact[i] = obj_idx;
    wall_strain[i]  = strain;

    // output collision energy
    if (is_first_collision && calc_energy) {
        if (step % colli_output_freq == 0)
        {
            double colllision_en = 0.5 * mass[species_cr[i]] * (rel_vel_n.x * rel_vel_n.x + rel_vel_n.y * rel_vel_n.y + rel_vel_n.z * rel_vel_n.z);
            double colllision_et = 0.5 * mass[species_cr[i]] * (rel_vel_t.x * rel_vel_t.x + rel_vel_t.y * rel_vel_t.y + rel_vel_t.z * rel_vel_t.z);
            fprintf(fp_energy_obj, "%8d,%8d,%8d,% .4e,%.4e,%.4e,%.4e,%.4e\n",
                step, i, obj_idx, wall_pos.x, wall_pos.y, wall_pos.z, colllision_en, colllision_et);
        }
    }
}

void DEM::calculation()
{
    // reset force and torque
    for (int i = 0; i < num_particle; ++i) {
        force[i].initialize();

        torque[i].initialize();
    }

    // particle - particle collision
    build_near_particle_list();
   
    for (int i = 0; i < num_particle; ++i)
    {
        int num_contact = 0;
        for (int jj = species_cr[i]; jj < num_species; ++jj)
        {
            const int check_pid = celldata[jj].cell_number[i] / (celldata[jj].num_cell_x * celldata[jj].num_cell_y);

            if (check_pid == 1) {
                // periodic boundary min : collision detection
                const int pid_x = celldata[jj].cell_number[i] % celldata[jj].num_cell_x;
                const int pid_y = (celldata[jj].cell_number[i] / celldata[jj].num_cell_x) % celldata[jj].num_cell_y;
                const int c_idx = pid_x + celldata[jj].num_cell_x * pid_y + celldata[jj].num_cell_x * celldata[jj].num_cell_y * (celldata[jj].num_cell_z - 2);

                for (int y = -1; y <= 1; ++y) {
                    const int pid = c_idx + celldata[jj].num_cell_x * y;
                    const int start = celldata[jj].count_per_cell[pid - 1];
                    const int end = celldata[jj].count_per_cell[pid + 2];

                    for (int k = start; k < end; ++k) {
                        const int j = celldata[jj].sort_list[k];

                        if (species_cr[j] != jj)
                        {
                            continue;
                        }

                        if (i < j) {
                            Vector3d norm;
                            norm.x = pos[i].x - pos[j].x;
                            norm.y = pos[i].y - pos[j].y;
                            norm.z = pos[i].z - (pos[j].z - periodic_dist);
                            double dist = norm.norm();
                            double distance_bet_surface = dist - radius[species_cr[i]] - radius[species_cr[j]];

                            if (calc_VanDerWaals_pp)
                            {
                                if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
                                {
                                    norm.normalize();
                                    calculate_VanDerWaals_pp(i, j, distance_bet_surface, norm);
                                }
                            }

                            if (distance_bet_surface < 0) {
                                norm.normalize();
                                calculate_ptcl_collision(i, j, num_contact, dist, norm);
                            }
                        }
                    }
                }

                for (int z = 0; z <= 1; ++z) {
                    for (int y = -1; y <= 1; ++y) {
                        const int pid = celldata[jj].cell_number[i] + celldata[jj].num_cell_x * y + celldata[jj].num_cell_x * celldata[jj].num_cell_y * z;
                        const int start = celldata[jj].count_per_cell[pid - 1];
                        const int end = celldata[jj].count_per_cell[pid + 2];

                        for (int k = start; k < end; ++k) {
                            const int j = celldata[jj].sort_list[k];

                            if (species_cr[j] != jj)
                            {
                                continue;
                            }

                            if (i < j) {
                                Vector3d norm = pos[i] - pos[j];
                                double   dist = norm.norm();
                                double distance_bet_surface = dist - radius[species_cr[i]] - radius[species_cr[j]];

                                if (calc_VanDerWaals_pp)
                                {
                                    if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
                                    {
                                        norm.normalize();
                                        calculate_VanDerWaals_pp(i, j, distance_bet_surface, norm);
                                    }
                                }

                                if (distance_bet_surface < 0) {
                                    norm.normalize();
                                    calculate_ptcl_collision(i, j, num_contact, dist, norm);
                                }
                            }
                        }
                    }
                }
            }
            else if (check_pid == (celldata[jj].num_cell_z - 2)) {
                // periodic boundary max : collision detection
                for (int z = -1; z <= 0; ++z) {
                    for (int y = -1; y <= 1; ++y) {
                        const int pid = celldata[jj].cell_number[i] + celldata[jj].num_cell_x * y + celldata[jj].num_cell_x * celldata[jj].num_cell_y * z;
                        const int start = celldata[jj].count_per_cell[pid - 1];
                        const int end = celldata[jj].count_per_cell[pid + 2];

                        for (int k = start; k < end; ++k) {
                            const int j = celldata[jj].sort_list[k];

                            if (species_cr[j] != jj)
                            {
                                continue;
                            }

                            if (i < j) {
                                Vector3d norm = pos[i] - pos[j];
                                double   dist = norm.norm();
                                double distance_bet_surface = dist - radius[species_cr[i]] - radius[species_cr[j]];

                                if (calc_VanDerWaals_pp)
                                {
                                    if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
                                    {
                                        norm.normalize();
                                        calculate_VanDerWaals_pp(i, j, distance_bet_surface, norm);
                                    }
                                }

                                if (distance_bet_surface < 0) {
                                    norm.normalize();
                                    calculate_ptcl_collision(i, j, num_contact, dist, norm);
                                }
                            }
                        }
                    }
                }

                const int pid_x = celldata[jj].cell_number[i] % celldata[jj].num_cell_x;
                const int pid_y = (celldata[jj].cell_number[i] / celldata[jj].num_cell_x) % celldata[jj].num_cell_y;
                const int c_idx = pid_x + celldata[jj].num_cell_x * pid_y + celldata[jj].num_cell_x * celldata[jj].num_cell_y;     // pid_z = 1
                for (int y = -1; y <= 1; ++y) {
                    const int pid = c_idx + celldata[jj].num_cell_x * y;
                    const int start = celldata[jj].count_per_cell[pid - 1];
                    const int end = celldata[jj].count_per_cell[pid + 2];

                    for (int k = start; k < end; ++k) {
                        const int j = celldata[jj].sort_list[k];

                        if (species_cr[j] != jj)
                        {
                            continue;
                        }

                        if (i < j) {
                            Vector3d norm;
                            norm.x = pos[i].x - pos[j].x;
                            norm.y = pos[i].y - pos[j].y;
                            norm.z = pos[i].z - (pos[j].z + periodic_dist);
                            double dist = norm.norm();
                            double distance_bet_surface = dist - radius[species_cr[i]] - radius[species_cr[j]];

                            if (calc_VanDerWaals_pp)
                            {
                                if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
                                {
                                    norm.normalize();
                                    calculate_VanDerWaals_pp(i, j, distance_bet_surface, norm);
                                }
                            }

                            if (distance_bet_surface < 0) {
                                norm.normalize();
                                calculate_ptcl_collision(i, j, num_contact, dist, norm);
                            }
                        }
                    }
                }
            }
            else {
                // normal collision detection
                for (int z = -1; z <= 1; ++z) {
                    for (int y = -1; y <= 1; ++y) {
                        const int pid = celldata[jj].cell_number[i] + celldata[jj].num_cell_x * y + celldata[jj].num_cell_x * celldata[jj].num_cell_y * z;
                        const int start = celldata[jj].count_per_cell[pid - 1];
                        const int end = celldata[jj].count_per_cell[pid + 2];

                        for (int k = start; k < end; ++k) {
                            const int j = celldata[jj].sort_list[k];

                            if (species_cr[j] != jj)
                            {
                                continue;
                            }

                            if (i < j) {
                                Vector3d norm = pos[i] - pos[j];
                                double   dist = norm.norm();
                                double distance_bet_surface = dist - radius[species_cr[i]] - radius[species_cr[j]];

                                if (calc_VanDerWaals_pp)
                                {
                                    if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
                                    {
                                        norm.normalize();
                                        calculate_VanDerWaals_pp(i, j, distance_bet_surface, norm);
                                    }
                                }

                                if (distance_bet_surface < 0) {
                                    norm.normalize();
                                    calculate_ptcl_collision(i, j, num_contact, dist, norm);
                                }
                            }
                        }
                    }
                }
            }
        }
 
        // update contact list
        ptcl_contact_num[i] = num_contact;
        for (int n = 0; n < 12; ++n) {  // num_possible_contact = 12
            ptcl_contact_list_prev[i][n] = ptcl_contact_list[i][n];
            ptcl_contact_strain_prev[i][n] = ptcl_contact_strain[i][n];

            ptcl_contact_list[i][n] = -1;
            ptcl_contact_strain[i][n].initialize();
        } 
    }
    
    //print contact number
    if (mes_contact_num || k_energy)
    {
        print_output();
    }
 
    // particle - wall collision
    for (int i = 0; i < num_particle; ++i) {
        int      obj_idx = -1;
        Vector3d wall_norm;
        
        double dist = wall->get_obj_distance(pos[i], radius[species_cr[i]], obj_idx, wall_norm);

        if (calc_VanDerWaals_pw)
        {
            double distance_bet_surface = dist - radius[species_cr[i]];
          
            if (cutoff_dist_min < distance_bet_surface && distance_bet_surface < cutoff_dist_max)
            {
                calculate_VanDerWaals_pw(i, dist, wall_norm);
            }
        }

        if (dist < radius[species_cr[i]]) {
            calculate_wall_collision(i, obj_idx, dist, wall_norm);
        }
    }

    //calculate kinetic energy
    if (k_energy)
    {
        calculate_ptcl_kinetic_energy();
    }

    // update
    for (int i = 0; i < num_particle; ++i) {
        vel[i] += (force[i] / mass[species_cr[i]] + gravity) * dt;

        pos[i] += vel[i] * dt;

        // z-direction periodic boundary
        if (periodic_max <= pos[i].z) {
            pos[i].z -= periodic_dist;
        }
        else if (pos[i].z < periodic_min) {
            pos[i].z += periodic_dist;
        }

        omega[i] += (torque[i] / moment[species_cr[i]]) * dt;

    }
 
    wall->update();
    step += 1;
}

