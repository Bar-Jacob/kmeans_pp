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

static PyObject *fit(PyObject *self, PyObject *args);
static int kmeans(int k, int max_iter, Py_ssize_t dimension_p, Py_ssize_t num_of_points_p, PyObject* centroids_locations, PyObject* data_points_p);
double Euclidian_Distance(double *vector1, double *vector2, int dimension);
void finding_cluster(double *vector, Cluster *clusters, int k, int dimension);
int update_mean(Cluster *clusters, int same_average, int k, int dimension);
void update_sum_of_elements_in_cluster(double *vector, int loc, Cluster *clusters, int dimension);
void free_memory(Cluster *clusters, double **data_points, int k, int num_of_points);

/*static int hello(int n, double z)
{
    printf("hi %d %f", n, z);
    return z;
}*/

static PyObject *fit(PyObject *self, PyObject *args)
{
    int k;
    int max_iter;
    Py_ssize_t dimension_p;
    Py_ssize_t num_of_points_p;
    Py_ssize_t n;
    PyObject *centroids_locations;
    PyObject *data_points_p;

    if (!PyArg_ParseTuple(args, "iiO!O!O!", &k, &max_iter, &dimension_p, &centroids_locations, &data_points_p))
    {
        return NULL;
    }
    if(!PYList_Check(centroids_locations) || !PYList_Check(data_points_p)){
        return NULL;
    }
    n = PyList_Size(data_points_p);
    num_of_points_p = n/dimension_p;
    kmeans(k, max_iter, dimension_p, num_of_points_p, centroids_locations, data_points_p);
    return PY_BuildValue(1);
}

static PyMethodDef kmeansMethods[] = {
    {"kmeans",
     (PyCFunction)fit,
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

static int kmeans(int k, int max_iter, Py_ssize_t dimension_p, Py_ssize_t num_of_points_p, PyObject* centroids_locations, PyObject* data_points_p)
{
    Cluster *clusters;
    double **data_points;
    int same_average = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    Py_ssize_t p = 0;
    Py_ssize_t q = 0;
    Py_ssize_t cnt_p = 0;
    int num_of_points = 0;
    int dimension = 0;


    num_of_points = (int)num_of_points_p;
    dimension = (int)dimension_p;

    data_points = (double **)calloc(num_of_points_p, sizeof(*data_points));
    assert(data_points != NULL);

    for (p = 0; p < num_of_points_p; p++)
    {
        for(q = 0; q < dimension_p; q++){
            data_points[p] = (double *)calloc(dimension_p, sizeof(*data_points[p]));
            assert(data_points[p] != NULL);
            data_points[p][q] = (double)PyList_GetItem(data_points_p, cnt_p);
            cnt_p++;
        }
    }

    /*
    Initializing k clusters
    */
    cnt_p = 0;
    clusters = (Cluster *)calloc(k, sizeof(struct Cluster));
    for (i = 0; i < k; i++)
    {
        clusters[i].centroid = (double *)calloc(dimension_p, sizeof(double));
        assert(clusters[i].centroid != NULL);

        memcpy(clusters[i].centroid, data_points[PyList_GetItem(centroids_locations,cnt_p)], sizeof(double) * dimension); /*will be equal to the i'th vector
        */
        clusters[i].num_of_points = 0;
        clusters[i].sum_of_points = (double *)calloc(dimension_p, sizeof(double));
        assert(clusters[i].sum_of_points != NULL);
        cnt_p++;
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
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            if (j == dimension - 1)
            {
                printf("%.4f", clusters[i].centroid[j]);
            }
            else
            {
                printf("%.4f,", clusters[i].centroid[j]);
            }
        }
        printf("\n");
    }
    free_memory(clusters, data_points, k, num_of_points);

    return 0;
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
