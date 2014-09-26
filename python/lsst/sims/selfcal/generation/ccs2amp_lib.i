%module ccs2amp_lib
%include "cpointer.i"
%include "carrays.i"
%include "cstring.i"

%pointer_functions(int, intp);
%array_functions(int, intA);
%array_functions(float, floatA);

%{
#define SWIG_FILE_WITH_INIT
#include "ampcoordtran.h"
%}

typedef struct {
int   raft_ix[2];
int   sens_ix[2];
int   device_ix[2];
float xo[2];
float xc[2][2];
float po[2];
float pc[2][2];
int   addr_lim[2][2];
float addr_elim[2][2];
} amp_coordtran;

int ccs2amp(double xFp, double yFp, amp_coordtran *ampDesc, int *iPix, int *jPix);

%cstring_bounded_output(char *detName, 16);
int ccs2ampbrief(double xFp, double yFp, char *detName, int *iPix, int *jPix);