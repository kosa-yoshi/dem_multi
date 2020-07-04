#ifndef DEM_CFD_PROJECT_WALL_BOUNDARY_H
#define DEM_CFD_PROJECT_WALL_BOUNDARY_H

#define _CRT_SECURE_NO_WARNINGS

#include <float.h>
#include "utility.h"

class Line
{
public:
    int index[2];
    double length;
    Vector3d direction;
};

class Triangle
{
public:
    int index[3];
    Vector3d vec_ab, vec_ac, norm;
};

class Axis
{
public:
    double rps;
    double rotation_matrix[3][3];
    Vector3d norm;
    Vector3d pos;

    Axis();
};

class WallObject
{
public:
    int num_point;
    int num_line;
    int num_mesh;
    int num_axis;

    WallObject();
    ~WallObject();

    bool load_obj_parameter(const char* file_name);
    void allocate_rotation_axis();
    void initialize_axis_parameter(const int id, const double rpm, const Vector3d norm, const Vector3d pos);
    void calculate_rotation_matrix(const double dt);

    void calculate_mesh_vector();
    void calculate_line_vector();

    double get_nearest_distance(const Vector3d& ptcl_pos, const double radius, Vector3d& v_norm);
    void   get_wall_velocity(const Vector3d& pos, Vector3d& vel);

    void get_min_point_position(Vector3d& min_pos);
    void get_max_point_position(Vector3d& max_pos);

    void update();
    void write_vtk_file(const int output_idx, const int obj_id);

private:
    Vector3d* point;
    Line* line;
    Triangle* mesh;
    Axis* axis;
};

class WallManager
{
public:
    WallObject* obj;

    WallManager();
    ~WallManager();

    void initialize();
    void set_object_param(const double dt);
    void write_vtk_file(const int output_idx);
    void update();

    double get_obj_distance(const Vector3d& ptcl_pos, const double radius, int& obj_idx, Vector3d& w_norm);

    void get_min_point_position(Vector3d& min_pos);
    void get_max_point_position(Vector3d& max_pos);

    int  get_num_object();


private:
    int num_object;
};

#endif
