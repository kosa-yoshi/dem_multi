#include "utility.h"

// -- class Util functions definition -- //
FILE* Util::file_check_open(const char *name, const char *mode)
{
    FILE *fp;
    fp  = fopen(name, mode);
    if (fp == NULL) {
        Log::write("[ERROR] *** file open failed [%s] ***\n", name);
        exit(1);
    }

    return fp;
}

void Util::get_int_data(FILE *fp, const char *name, int& data)
{
    char buf[128], tag[64];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %d", tag, &data);
    Log::write("%s:: %d\n", name, data);
}

void Util::get_double_data(FILE *fp, const char *name, double& data)
{
    char buf[128], tag[64];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %lf", tag, &data);
    Log::write("%s:: % .4e\n", name, data);
}

void Util::get_char_data(FILE *fp, const char *name, char *string)
{
    char buf[128], tag[64];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %s", tag, string);
    Log::write("%s:: %s\n", name, string);
}

void Util::get_bool_data(FILE *fp, const char *name, bool& flag)
{
    char buf[128], tag[64], val[16];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %s", tag, val);
    if ((strstr(val, "yes")) != NULL) {
        flag = true;
        Log::write("%s : yes\n", name);
    } else {
        flag = false;
        Log::write("%s : no\n", name);
    }
}

void Util::get_vec_data(FILE *fp, const char *name, Vector3d& vec)
{
    char buf[128], tag[64];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %lf %lf %lf", tag, &vec.x, &vec.y, &vec.z);
    Log::write("%s", buf);
}

void Util::get_vec_data(FILE *fp, Vector3d& vec)
{
    char buf[128];
    fgets(buf, sizeof(buf), fp);
    sscanf(buf, "%lf %lf %lf", &vec.x, &vec.y, &vec.z);
}

void Util::get_int3_data(FILE *fp, const char *name, int& i, int& j, int& k)
{
    char buf[128], tag[64];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    sscanf(buf, "%s %d %d %d", tag, &i, &j, &k);
    Log::write("%s", buf);
}

void Util::get_int3_data(FILE *fp, int& i, int& j, int& k)
{
    char buf[128];
    fgets(buf, sizeof(buf), fp);
    sscanf(buf, "%d %d %d", &i, &j, &k);
}

void Util::get_int2_data(FILE *fp, int& i, int& j)
{
    char buf[128];
    fgets(buf, sizeof(buf), fp);
    sscanf(buf, "%d %d", &i, &j);
}

void Util::read_line_data(FILE *fp, const char *name)
{
    char buf[128];
    fgets(buf, sizeof(buf), fp);
    validate_read_data(buf, name);
    Log::write("%s", buf);
}

void Util::validate_read_data(const char *buf, const char *name)
{
    if ((strstr(buf, name)) == NULL) {
        Log::write("[error] unexpected data in file\n>> %s", buf);
        Log::write(">> [%s] should be expected\n", name);
        Log::close();
        exit(1);
    }
}

bool Util::check_line_data(FILE *fp, const char *name)
{
    char buf[128];
    fgets(buf, sizeof(buf), fp);
    if ((strstr(buf, name)) == NULL) {
        Log::write("[%s] data dose not exist in file\n", name);
        return false;
    } else {
        Log::write("[%s] data exist in file\n", name);
        return true;
    }
        return false;
}

void Util::show_vector_data(const Vector3d& vec)
{
    Log::write("% .4e, % .4e, % .4e\n", vec.x, vec.y, vec.z);
}

void Util::show_vector_data(const char *name, const Vector3d& vec)
{
    Log::write("[%s] : % .4e, % .4e, % .4e\n", name, vec.x, vec.y, vec.z);
}

void Util::allocate_vector_array_3d(Vector3d*** &ptr, int nx, int ny, int nz)
{
    ptr       = new Vector3d** [nx];
    ptr[0]    = new Vector3d*  [nx * ny];
    ptr[0][0] = new Vector3d   [nx * ny * nz];

    for (int i = 0; i < nx; ++i) {
        ptr[i] = ptr[0] + i * ny;
        for (int j = 0; j < ny; ++j) {
            ptr[i][j] = ptr[0][0] + i * ny * nz + j * nz;
        }
    }
}

void Util::allocate_vector_array_2d(Vector3d** &ptr, int nx, int ny)
{
    ptr    = new Vector3d* [nx];
    ptr[0] = new Vector3d  [nx * ny];

    for (int i = 0; i < nx; ++i) {
        ptr[i] = ptr[0] + i * ny;
    }
}

void Util::allocate_double_array_3d(double*** &ptr, int nx, int ny, int nz)
{
    ptr       = new double** [nx];
    ptr[0]    = new double*  [nx * ny];
    ptr[0][0] = new double   [nx * ny * nz];

    for (int i = 0; i < nx; ++i) {
        ptr[i] = ptr[0] + i * ny;
        for (int j = 0; j < ny; ++j) {
            ptr[i][j] = ptr[0][0] + i * ny * nz + j * nz;
        }
    }
}

void Util::allocate_double_array_2d(double** &ptr, int nx, int ny)
{
    ptr    = new double* [nx];
    ptr[0] = new double  [nx * ny];

    for (int i = 0; i < nx; ++i) {
        ptr[i] = ptr[0] + i * ny;
    }
}

void Util::allocate_int_array_2d(int** &ptr, int nx, int ny)
{
    ptr    = new int* [nx];
    ptr[0] = new int  [nx * ny];

    for (int i = 0; i < nx; ++i) {
        ptr[i] = ptr[0] + i * ny;
    }
}

void Util::free_vector_array_3d(Vector3d*** &ptr)
{
    delete [] ptr[0][0];
    delete [] ptr[0];
    delete [] ptr;
}

void Util::free_vector_array_2d(Vector3d** &ptr)
{
    delete [] ptr[0];
    delete [] ptr;
}

void Util::free_double_array_3d(double*** &ptr)
{
    delete [] ptr[0][0];
    delete [] ptr[0];
    delete [] ptr;
}

void Util::free_double_array_2d(double** &ptr)
{
    delete [] ptr[0];
    delete [] ptr;
}

void Util::free_int_array_2d(int** &ptr)
{
    delete [] ptr[0];
    delete [] ptr;
}

void Util::make_transposed_matrix(double inp[][3], double out[][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i][j] = inp[j][i];
        }
    }
}

void Util::multi_matrix_matrix(double left[][3], double right[][3], double out[][3])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                out[i][j] += left[i][k] * right[k][j];
            }
        }
    }
}

void Util::multi_matrix_vector(double matrix[][3], Vector3d& vec_inp, Vector3d& vec_out)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vec_out[i] += matrix[i][j] * vec_inp[j];
        }
    }
}

// -- class Log definition -- //
FILE* Log::fp;
bool  Log::log;

void Log::open(const char* name, const char* title)
{
    fp = Util::file_check_open(name, "w");
    log = true;

    time_t t;
    time(&t);
    write("%s : %s\n", title, ctime(&t));
}

void Log::close()
{
    fclose(fp);
    log = false;
}

int Log::write(const char* format, ...)
{
    va_list list;
    va_start(list, format);

    char buf[256];
    int  num = vsprintf(buf, format, list);
    fputs(buf, stdout);
    fflush(stdout);

    if (log) {
        fputs(buf, fp);
        fflush(fp);
    }

    va_end(list);
    return num;
}

int Log::write_file(const char* format, ...)
{
    va_list list;
    va_start(list, format);

    char buf[256];
    int  num = vsprintf(buf, format, list);

    if (log) {
        fputs(buf, fp);
        fflush(fp);
    }

    va_end(list);
    return num;
}

