#include "wall_boundary.h"

Axis::Axis()
{
    rps = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rotation_matrix[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

WallObject::WallObject()
{
    point = nullptr;
    line = nullptr;
    mesh = nullptr;
    axis = nullptr;

    num_point = 0;
    num_line = 0;
    num_mesh = 0;
    num_axis = 0;
}

WallObject::~WallObject()
{
    delete[] point;
    delete[] line;
    delete[] mesh;
    delete[] axis;
}

bool WallObject::load_obj_parameter(const char* file_name)
{
    FILE* fp;
    fp = Util::file_check_open(file_name, "r");
    Util::get_int_data(fp, "Number_of_Point", num_point);
    try {
        point = new Vector3d[num_point];
    }
    catch (std::bad_alloc) {
        Log::write("[ERROR] *** memory exhausted (point) ***\n");
        return false;
    }

    for (int i = 0; i < num_point; ++i) {
        Util::get_vec_data(fp, point[i]);
    }

    Util::get_int_data(fp, "Number_of_Mesh", num_mesh);
    try {
        mesh = new Triangle[num_mesh];
    }
    catch (std::bad_alloc) {
        Log::write("[ERROR] *** memory exhausted (mesh) ***\n");
        return false;
    }

    for (int i = 0; i < num_mesh; ++i) {
        Util::get_int3_data(fp, mesh[i].index[0], mesh[i].index[1], mesh[i].index[2]);
    }

    Util::get_int_data(fp, "Number_of_Line", num_line);
    try {
        line = new Line[num_line];
    }
    catch (std::bad_alloc) {
        Log::write("[ERROR] *** memory exhausted (line) ***\n");
        return false;
    }

    for (int i = 0; i < num_line; ++i) {
        Util::get_int2_data(fp, line[i].index[0], line[i].index[1]);
    }

    fclose(fp);
    return true;
}

void WallObject::allocate_rotation_axis()
{
    axis = new Axis[num_axis];
}

void WallObject::initialize_axis_parameter(const int id, const double rpm, const Vector3d norm, const Vector3d pos)
{
    axis[id].rps = rpm / 60.0;

    axis[id].norm = norm;
    axis[id].norm.normalize();

    axis[id].pos = pos;
}

void WallObject::calculate_rotation_matrix(const double dt)
{
    for (int i = 0; i < num_axis; ++i) {
        double theta = 2.0 * 3.14 * axis[i].rps * dt;
        double cos_dt = cos(theta);
        double sin_dt = sin(theta);

        axis[i].rotation_matrix[0][0] = axis[i].norm.x * axis[i].norm.x * (1.0 - cos_dt) + cos_dt;
        axis[i].rotation_matrix[0][1] = axis[i].norm.x * axis[i].norm.y * (1.0 - cos_dt) + axis[i].norm.z * sin_dt;
        axis[i].rotation_matrix[0][2] = axis[i].norm.x * axis[i].norm.z * (1.0 - cos_dt) - axis[i].norm.y * sin_dt;

        axis[i].rotation_matrix[1][0] = axis[i].norm.y * axis[i].norm.x * (1.0 - cos_dt) - axis[i].norm.z * sin_dt;
        axis[i].rotation_matrix[1][1] = axis[i].norm.y * axis[i].norm.y * (1.0 - cos_dt) + cos_dt;
        axis[i].rotation_matrix[1][2] = axis[i].norm.y * axis[i].norm.z * (1.0 - cos_dt) + axis[i].norm.x * sin_dt;

        axis[i].rotation_matrix[2][0] = axis[i].norm.z * axis[i].norm.x * (1.0 - cos_dt) + axis[i].norm.y * sin_dt;
        axis[i].rotation_matrix[2][1] = axis[i].norm.z * axis[i].norm.y * (1.0 - cos_dt) - axis[i].norm.x * sin_dt;
        axis[i].rotation_matrix[2][2] = axis[i].norm.z * axis[i].norm.z * (1.0 - cos_dt) + cos_dt;
    }
}

void WallObject::calculate_mesh_vector()
{
    for (int i = 0; i < num_mesh; ++i) {
        mesh[i].vec_ab = point[mesh[i].index[1]] - point[mesh[i].index[0]];
        mesh[i].vec_ac = point[mesh[i].index[2]] - point[mesh[i].index[0]];
        mesh[i].norm = mesh[i].vec_ac % mesh[i].vec_ab;
        mesh[i].norm.normalize();
    }
}

void WallObject::calculate_line_vector()
{
    for (int i = 0; i < num_line; ++i) {
        line[i].direction = point[line[i].index[1]] - point[line[i].index[0]];
        line[i].length = line[i].direction.norm();
        line[i].direction.normalize();
    }
}

double WallObject::get_nearest_distance(const Vector3d& ptcl_pos, const double radius, Vector3d& v_norm)
{
    double min_dist = DBL_MAX;
    const double eps = DBL_EPSILON * 10.0;
    for (int i = 0; i < num_mesh; ++i) {
        Vector3d vec = ptcl_pos - point[mesh[i].index[0]];
        double dist = fabs(vec * mesh[i].norm);

        if (dist > 2 * radius) {
            continue;
        }

        Vector3d p = vec + dist * mesh[i].norm;

        double pb = p * mesh[i].vec_ab;
        double pc = p * mesh[i].vec_ac;

        double bc = mesh[i].vec_ab * mesh[i].vec_ac;
        double bb = mesh[i].vec_ab * mesh[i].vec_ab;
        double cc = mesh[i].vec_ac * mesh[i].vec_ac;

        double div = bc * bc - bb * cc;
        double t = (pb * bc - pc * bb) / div;

        if (t <= eps || ((1.0 - eps) <= t)) {
            continue;
        }

        double s = (pc * bc - pb * cc) / div;
        if (s <= eps || ((1.0 - eps) <= (s + t))) {
            continue;
        }

        if (min_dist > dist) {
            min_dist = dist;
            v_norm = mesh[i].norm;
        }
    }

    for (int i = 0; i < num_line; ++i) {
        Vector3d vec = ptcl_pos - point[line[i].index[0]];
        double   dot = vec * line[i].direction;
        if (dot <= eps || line[i].length <= (dot + eps)) {
            continue;
        }

        Vector3d norm = dot * line[i].direction - vec;
        double   dist = norm.norm();
        norm.normalize();
        norm *= -1.0;

        if (dist > 2 * radius) {
            continue;
        }

        if (dist < min_dist) {
            min_dist = dist;
            v_norm = norm;
        }
    }

    return min_dist;
}

void WallObject::get_wall_velocity(const Vector3d& pos, Vector3d& vel)
{
    for (int i = 0; i < num_axis; ++i) {
        vel += 2.0 * 3.14 * axis[i].rps * ((pos - axis[i].pos) % axis[i].norm);
    }
}

void WallObject::get_min_point_position(Vector3d& min_pos)
{
    for (int i = 0; i < num_point; ++i) {
        min_pos.x = (point[i].x < min_pos.x) ? point[i].x : min_pos.x;
        min_pos.y = (point[i].y < min_pos.y) ? point[i].y : min_pos.y;
        min_pos.z = (point[i].z < min_pos.z) ? point[i].z : min_pos.z;
    }
}

void WallObject::get_max_point_position(Vector3d& max_pos)
{
    for (int i = 0; i < num_point; ++i) {
        max_pos.x = (point[i].x > max_pos.x) ? point[i].x : max_pos.x;
        max_pos.y = (point[i].y > max_pos.y) ? point[i].y : max_pos.y;
        max_pos.z = (point[i].z > max_pos.z) ? point[i].z : max_pos.z;
    }
}

void WallObject::update()
{
    for (int i = 0; i < num_axis; ++i) {
        for (int j = 0; j < num_point; ++j) {
            Vector3d pre_pos = point[j] - axis[i].pos;
            Vector3d new_pos;

            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    new_pos[m] += axis[i].rotation_matrix[m][n] * pre_pos[n];
                }
            }
            point[j] = new_pos + axis[i].pos;
        }
    }

    calculate_mesh_vector();

    calculate_line_vector();
}

void WallObject::write_vtk_file(const int output_idx, const int obj_id)
{
    char name[64];
    sprintf(name, "./output/wall_obj_%02d_%04d.vtk", obj_id, output_idx);

    FILE* fp;
    fp = Util::file_check_open(name, "w");

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "wall_obj_data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");

    fprintf(fp, "POINTS %d double\n", num_point);
    for (int i = 0; i < num_point; ++i) {
        fprintf(fp, "% .4e % .4e % .4e\n", point[i].x, point[i].y, point[i].z);
    }

    fprintf(fp, "POLYGONS %d %d\n", num_mesh, num_mesh * 4);
    for (int i = 0; i < num_mesh; ++i) {
        fprintf(fp, "3 %d %d %d\n", mesh[i].index[0], mesh[i].index[1], mesh[i].index[2]);
    }
    fclose(fp);
}

WallManager::WallManager()
{
    obj = nullptr;
}

WallManager::~WallManager()
{
    delete[] obj;
}

void WallManager::initialize()
{
    FILE* fp_obj;
    fp_obj = Util::file_check_open("./input/config_obj.txt", "r");

    Util::get_int_data(fp_obj, "num_object", num_object);
    obj = new WallObject[num_object];

    for (int i = 0; i < num_object; ++i) {
        int  object_index = 0;
        char file_name[128], object_file[256];

        Util::get_int_data(fp_obj, "object_index", object_index);
        Util::get_char_data(fp_obj, "file_name", file_name);

        sprintf(object_file, "./input/%s", file_name);
        if (!obj[i].load_obj_parameter(object_file)) {
            Log::write("[ERROR] *** load obj text file [%s] ***\n", file_name);
            exit(1);
        }

        bool rotation = false;
        Util::get_bool_data(fp_obj, "rotation", rotation);
        if (rotation) {
            Util::get_int_data(fp_obj, "num_axis", obj[i].num_axis);
            if (obj[i].num_axis <= 0) {
                Log::write("num_axis is invalied value [%d]\n", obj[i].num_axis);
                exit(2);
            }

            obj[i].allocate_rotation_axis();

            for (int j = 0; j < obj[i].num_axis; ++j) {
                double rpm = 0.0;
                Vector3d axis_norm, axis_pos;
                Util::get_double_data(fp_obj, "rot_speed", rpm);
                Util::get_vec_data(fp_obj, "axis_norm", axis_norm);
                Util::get_vec_data(fp_obj, "axis_pos", axis_pos);

                obj[i].initialize_axis_parameter(j, rpm, axis_norm, axis_pos);
            }
        }

        while (true) {
            char str[256];
            fgets(str, sizeof(str), fp_obj);
            if ((strstr(str, "end")) != NULL) {
                break;
            }
        }
    }

    fclose(fp_obj);
}

void WallManager::set_object_param(const double dt)
{
    for (int i = 0; i < num_object; ++i) {
        obj[i].calculate_mesh_vector();
        obj[i].calculate_line_vector();

        if (obj[i].num_axis > 0) {
            obj[i].calculate_rotation_matrix(dt);
        }
    }
}

void WallManager::write_vtk_file(const int output_idx)
{
    for (int i = 0; i < num_object; ++i) {
        obj[i].write_vtk_file(output_idx, i);
    }
}

void WallManager::update()
{
    for (int i = 0; i < num_object; ++i) {
        if (obj[i].num_axis > 0) {
            obj[i].update();
        }
    }
}

double WallManager::get_obj_distance(const Vector3d& ptcl_pos, const double radius, int& obj_idx, Vector3d& w_norm)
{
    double obj_dist = DBL_MAX;
    for (int i = 0; i < num_object; ++i) {
        double distance = obj[i].get_nearest_distance(ptcl_pos, radius, w_norm);

        if (obj_dist > distance) {
            obj_dist = distance;
            obj_idx = i;
        }
    }

    return obj_dist;
}

void WallManager::get_min_point_position(Vector3d& min_pos)
{
    min_pos.set_value(DBL_MAX, DBL_MAX, DBL_MAX);

    for (int i = 0; i < num_object; ++i) {
        obj[i].get_min_point_position(min_pos);
    }
}

void WallManager::get_max_point_position(Vector3d& max_pos)
{
    max_pos.set_value(-DBL_MAX, -DBL_MAX, -DBL_MAX);

    for (int i = 0; i < num_object; ++i) {
        obj[i].get_max_point_position(max_pos);
    }
}

int WallManager::get_num_object()
{
    return num_object;
}
