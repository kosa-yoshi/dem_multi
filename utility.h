#define _CRT_SECURE_NO_WARNINGS

#ifndef DEM_CFD_PROJECT_UTILITY_H
#define DEM_CFD_PROJECT_UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include "vector3d.h"

class Util
{
    public:
        static FILE* file_check_open(const char *name, const char *mode);
        static void get_int_data(FILE *fp, const char *name, int& data);
        static void get_double_data(FILE *fp, const char *name, double& data);
        static void get_char_data(FILE *fp, const char *name, char *string);
        static void get_bool_data(FILE *fp, const char *name, bool& flag);
        static void get_vec_data(FILE *fp, const char *name, Vector3d& vec);
        static void get_vec_data(FILE *fp, Vector3d& vec);
        static void get_int3_data(FILE *fp, const char *name, int& i, int& j, int& k);
        static void get_int3_data(FILE *fp, int& i, int& j, int& k);
        static void get_int2_data(FILE *fp, int& i, int& j);
        static void read_line_data(FILE *fp, const char *name);
        static bool check_line_data(FILE *fp, const char *name);

        static void show_vector_data(const Vector3d& vec);
        static void show_vector_data(const char *name, const Vector3d& vec);

        static void allocate_vector_array_3d(Vector3d*** &ptr, int nx, int ny, int nz);
        static void allocate_vector_array_2d(Vector3d** &ptr, int nx, int ny);
        static void allocate_double_array_3d(double*** &ptr, int nx, int ny, int nz);
        static void allocate_double_array_2d(double** &ptr, int nx, int ny);
        static void allocate_int_array_2d(int** &ptr, int nx, int ny);

        static void free_vector_array_3d(Vector3d*** &ptr);
        static void free_vector_array_2d(Vector3d** &ptr);
        static void free_double_array_3d(double*** &ptr);
        static void free_double_array_2d(double** &ptr);
        static void free_int_array_2d(int** &ptr);

        static void make_transposed_matrix(double inp[][3], double out[][3]);
        static void multi_matrix_matrix(double left[][3], double right[][3], double out[][3]);
        static void multi_matrix_vector(double matrix[][3], Vector3d& vec_inp, Vector3d& vec_out);

    private:
        static void validate_read_data(const char *buf, const char *name);
};

class NonCopiable
{
    protected:
        NonCopiable() {}
        ~NonCopiable() {}

    private:
        void operator =(const NonCopiable& src);
        NonCopiable(const NonCopiable& src);
};

class Log : private NonCopiable
{
    public:
        Log() {log = false;}
        ~Log() {}

        static void open(const char* name, const char* title);
        static void close();
        static int write(const char* format, ...);
        static int write_file(const char* format, ...);

    private:
        static FILE *fp;
        static bool log;
};

#endif
