#include <Python.h>
/*
 * C extension module for python to count tetra nucleotides
 */

char int2res[4];
char int2rev[4];
int res2int['Z'+1];

void initialize(void)
{
  int2res[0] = 'A';
  int2res[1] = 'C';
  int2res[2] = 'T';
  int2res[3] = 'G';
  int2rev[0] = 2;
  int2rev[1] = 3;
  int2rev[2] = 0;
  int2rev[3] = 1;

  for (unsigned int i = 0; i <= 'Z'; ++i)
  {
    res2int[i]= -1;
  }
  for (unsigned int i = 0; i < 4; ++i)
  {
    res2int[(int)int2res[i]] = i;
  }
}

int *doTNCount(char *string){

    unsigned int length, valid=0;
    int *counts = calloc(256, sizeof(int));
    int nucl;
    
    length = strlen(string);
    uint8_t code = 0, revcode = 0;

    for(unsigned int i = 0; i < length; i++){
      
        nucl = res2int[(int)string[i]];
            
        if (nucl != -1)
        {
            valid++;
            code = (code << (uint8_t) 2) | (uint8_t)nucl;
            revcode = (revcode >> (uint8_t) 2) | 
                      ((uint8_t)int2rev[nucl] << (uint8_t) 6);

            if(valid >= 4){
                counts[code] += 1;
                counts[revcode] += 1;
            }
        }
        else
        {
            //reset values
            valid = 0;
            code = 0;
            revcode = 0;
        }
    }
    
    return counts;
}

static PyObject* count(PyObject* self, PyObject* args)
{
    char *str;
    int *counts;
    PyObject *count_list = PyList_New(256);

    if (!PyArg_ParseTuple(args, "s", &str))
        return NULL;

    counts = doTNCount(str);
    for (unsigned i = 0; i < 256; i++){
        PyList_SetItem(count_list, i, Py_BuildValue("i",counts[i]+1));
    }

    free(counts);
    return count_list;
}

// method mapping table
static PyMethodDef tnCounterMethods[] = {
    {"count", count, METH_VARARGS, "count tetranucleotides"},
    {NULL, NULL, 0, NULL}
};

// define module
static struct PyModuleDef tnCounterModule = {
    PyModuleDef_HEAD_INIT,
    "tnCounter",       // name of module
    NULL,              // module documentation, may be NULL
    -1,                // module keeps state in global variables
    tnCounterMethods
};

// init
PyMODINIT_FUNC PyInit_tnCounter(void)
{
    initialize();
    return PyModule_Create(&tnCounterModule);
}
