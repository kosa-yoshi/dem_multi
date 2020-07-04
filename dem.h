#ifndef DEM_CFD_PROJECT_DEM_H
#define DEM_CFD_PROJECT_DEM_H

#include "utility.h"
#include "wall_boundary.h"

class DEM
{
    public:
        int num_particle = 0;

        char axis[4];
        int  num_species;

        // particle property
        int*    num_particle_array;
        int*    species_cr;
        double* diameter;

        double*  radius;
        double*  volume;
        double*  density;
        double*  mass;
        double*  moment;
        double** dashpot;

        double** converted_diameter;
        double** converted_mass;

        double dt;
        double spring;
 
        double restitution;
        double friction, friction_wall;

        double ave_k_energy;
        double kinetic_energy;
        double kinetic_energy_pre;
        double min_k_energy;
        double max_k_energy;
        double hamaker_pp;
        double hamaker_pw;

        double cutoff_dist_min;
        double cutoff_dist_max;

        // periodic boundary
        double periodic_max;
        double periodic_min;
        double periodic_dist;

        //int**  contact_patern;

        int colli_output_freq;

        int count_00;
        int count_01;
        int count_11;

        //double cutoff_dist_pp_min;
        //double cutoff_dist_pp_max;
        //double cutoff_dist_pw_min;
        //double cutoff_dist_pw_max;

        // near particle search list
        struct Cell
        {
            double cell_size;
            int num_cell;
            int num_cell_x;
            int num_cell_y;
            int num_cell_z;

            Vector3d cell_pos;

            int* cell_number;
            int* count_per_cell;
            int* sort_list;

        } *celldata;

        Vector3d gravity;

        // particle data
        Vector3d* pos;
        Vector3d* vel;
        Vector3d* force;
        Vector3d* omega;
        Vector3d* torque;

        WallManager* wall;

        int*  ptcl_contact_num;
        int** ptcl_contact_list;
        int** ptcl_contact_list_prev;
        Vector3d** ptcl_contact_strain;
        Vector3d** ptcl_contact_strain_prev;

        int* wall_contact;
        Vector3d* wall_strain;


        DEM();
        ~DEM();

        void initialize();
        void finalization();

        void write_dem_output(const int output_id);
        void set_step(const int s);     // inline function
        void open_log_file(const bool b_energy, const int colli_frequency, const bool b_contct_num, 
            const bool b_kenergy, const bool b_van, const bool b_van_pw);
        void close_log_file();

        void calculation();
        void build_near_particle_list();
        void calculate_ptcl_collision(const int i, const int j, int& contact, const double dist, const Vector3d& norm);
        void calculate_wall_collision(const int i, const int obj_idx, const double dist, const Vector3d& wall_norm);


        void calculate_ptcl_kinetic_energy();
        void print_output();
        void calculate_VanDerWaals_pp(const int i, const int j, const double dist, const Vector3d& norm);
        void calculate_VanDerWaals_pw(const int i, const double dist, const Vector3d& wall_norm);
        

private:
    FILE* fp_over;

    FILE* fp_colli_ptcl;
    FILE* fp_colli_wall;

    FILE* fp_energy_ptcl, * fp_energy_obj;
    bool calc_energy;
    int step;

    FILE* fp_contact_num_ptcl;
    bool mes_contact_num;

    FILE* fp_kinetic_energy;
    bool k_energy;

    //FILE* fp_van_pp;
    bool calc_VanDerWaals_pp;

    //FILE* fp_van_pw;
    bool calc_VanDerWaals_pw;
    
};

inline void DEM::set_step(const int s)
{
    step = s;
}

#endif



        