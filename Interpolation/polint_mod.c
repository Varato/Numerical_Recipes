#include <Python.h>
#include "../nrlib/polint.h"

float interpolate(float x)
{
	int n=4;
	float xa[]={0, -1, 1, 2, 4};
	float ya[]={0, 1.25, 2, 3, 0};
	float y, dy;
	polint(xa, ya, n, x, &y, &dy);
	return y;
}

static PyObject * wrapped_interpolate(PyObject *self, PyObject *args)
{
	float _x;
	float y;
	if (!PyArg_ParseTuple(args, "d", &_x))
		return NULL;
	y=interpolate(_x);
	return Py_BuildValue("f", y);
}

static PyMethodDef polint_mod_methods[]={
	{"interpolate", wrapped_interpolate, METH_VARARGS, ""},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpolint_mod(void){
	(void) Py_InitModule("polint_mod", polint_mod_methods);
}
