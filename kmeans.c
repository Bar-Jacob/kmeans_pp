#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef struct Cluster
{
    int num_of_points;
    double *sum_of_points;
    double *centroid;
} Cluster;

static PyObject *fit_capi(PyObject *self, PyObject *args);
static PyObject *kmeans(int k, int max_iter, int dimension_p, int num_of_points_p, PyObject *centroids_locations, PyObject *data_points_p);
double Euclidian_Distance(double *vector1, double *vector2, int dimension);
void finding_cluster(double *vector, Cluster *clusters, int k, int dimension);
int update_mean(Cluster *clusters, int same_average, int k, int dimension);
void update_sum_of_elements_in_cluster(double *vector, int loc, Cluster *clusters, int dimension);
void free_memory(Cluster *clusters, double **data_points, int k, int num_of_points);
PyObject *cToPyObject(Cluster *clusters, int k, int dimension, double **data_points, int num_of_points);

int main(int argc, char *argv[])
{
    return 0;
}

static PyObject *kmeans(int k, int max_iter, int dimension_p, int num_of_points_p, PyObject *centroids_locations, PyObject *data_points_p)
{
    Cluster *clusters;
    double **data_points;
    int same_average = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    int num_of_points = num_of_points_p;
    int dimension = dimension_p;

    data_points = (double **)calloc(num_of_points, sizeof(*data_points));
    assert(data_points != NULL);

    for (i = 0; i < num_of_points; i++)
    {
        data_points[i] = (double *)calloc(dimension, sizeof(*data_points[i]));
        assert(data_points[i] != NULL);
        for (j = 0; j < dimension; j++)
        {
            data_points[i][j] = PyFloat_AsDouble(PyList_GetItem(data_points_p, cnt));
            cnt++;
        }
    }

    /*
    Initializing k clusters
    */
    cnt = 0;
    clusters = (Cluster *)calloc(k, sizeof(struct Cluster));
    assert(clusters != NULL);
    for (i = 0; i < k; i++)
    {
        clusters[i].centroid = (double *)calloc(dimension, sizeof(double));
        assert(clusters[i].centroid != NULL);

        memcpy(clusters[i].centroid, data_points[PyLong_AsLong(PyList_GetItem(centroids_locations, cnt))], sizeof(double) * dimension); /*will be equal to the i'th vector
        */
        clusters[i].num_of_points = 0;
        clusters[i].sum_of_points = (double *)calloc(dimension, sizeof(double));
        assert(clusters[i].sum_of_points != NULL);
        cnt++;
    }

    cnt = 0;
    while ((cnt < max_iter) && (!same_average))
    {
        same_average = 1;
        for (i = 0; i < num_of_points; i++)
        {
            finding_cluster(data_points[i], clusters, k, dimension);
        }

        same_average = update_mean(clusters, same_average, k, dimension);
        if (same_average == 1)
        {
            break;
        }

        for (i = 0; i < k; i++)
        {
            clusters[i].num_of_points = 0;
            for (j = 0; j < dimension; j++)
            {
                clusters[i].sum_of_points[j] = 0;
            }
        }
        cnt++;
    }

    return cToPyObject(clusters, k, dimension, data_points, num_of_points);
}

void finding_cluster(double *vector, Cluster *clusters, int k, int dimension)
{
    double min_distance = -1.0;
    int num_of_cluster = -1;
    double distance;
    int i = 0;

    for (i = 0; i < k; i++)
    {
        distance = Euclidian_Distance(vector, clusters[i].centroid, dimension);
        if ((distance < min_distance) || (min_distance < 0))
        {
            min_distance = distance;
            num_of_cluster = i;
        }
    }
    clusters[num_of_cluster].num_of_points++;
    update_sum_of_elements_in_cluster(vector, num_of_cluster, clusters, dimension);
}

void update_sum_of_elements_in_cluster(double *vector, int loc, Cluster *clusters, int dimension)
{
    int i = 0;
    for (i = 0; i < dimension; i++)
    {
        clusters[loc].sum_of_points[i] += vector[i];
    }
}

int update_mean(Cluster *clusters, int same_average, int k, int dimension)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            if ((clusters[i].sum_of_points[j] / clusters[i].num_of_points) !=
                clusters[i].centroid[j])
            {
                same_average = 0;
                clusters[i].centroid[j] = clusters[i].sum_of_points[j] / clusters[i].num_of_points;
            }
        }
    }
    return same_average;
}

double Euclidian_Distance(double *vector1, double *centroid, int dimension)
{
    double sum = 0.0;
    int xi = 0;
    for (xi = 0; xi < dimension; xi++)
    {
        sum += (vector1[xi] - centroid[xi]) * (vector1[xi] - centroid[xi]);
    }
    return sum;
}

PyObject *cToPyObject(Cluster *clusters, int k, int dimension, double **data_points, int num_of_points)
{
    PyObject *clusters_py;
    int i = 0;
    int j = 0;
    PyObject *value;

    clusters_py = PyList_New(k);
    for (i = 0; i < k; i++)
    {
        PyObject *curr_vector;
        curr_vector = PyList_New(dimension);
        for (j = 0; j < dimension; j++)
        {
            value = Py_BuildValue("d", clusters[i].centroid[j]);
            PyList_SetItem(curr_vector, j, value);
        }
        /*
        adding the PyObject centroid to the PyList clusters
        */
        PyList_SetItem(clusters_py, i, curr_vector);
    }
    free_memory(clusters, data_points, k, num_of_points);
    return clusters_py;
}

void free_memory(Cluster *clusters, double **data_points, int k, int num_of_points)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < num_of_points; i++)
    {
        free(data_points[i]);
    }
    free(data_points);

    for (j = 0; j < k; j++)
    {
        free(clusters[j].centroid);
        free(clusters[j].sum_of_points);
    }
    free(clusters);
}

static PyObject *fit_capi(PyObject *self, PyObject *args)
{
    int k;
    int max_iter;
    int dimension_p;
    int num_of_points_p;
    PyObject *centroids_locations;
    PyObject *data_points_p;

    if (!(PyArg_ParseTuple(args, "iiiiOO", &k, &max_iter, &dimension_p, &num_of_points_p, &centroids_locations, &data_points_p)))
    {
        return NULL;
    }
    if (!PyList_Check(centroids_locations) || !PyList_Check(data_points_p))
    {
        return NULL;
    }

    return Py_BuildValue("O", kmeans(k, max_iter, dimension_p, num_of_points_p, centroids_locations, data_points_p));
}

static PyMethodDef kmeansMethods[] = {
    {"fit",
     (PyCFunction)fit_capi,
     METH_VARARGS,
     PyDoc_STR("kmeans algorithem")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef =
    {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        kmeansMethods};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}