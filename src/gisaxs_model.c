#include <math.h>
//#include <R.h>

void ff_sphere(double *q, double *r, double *ff)
{
	double qr = *q * *r;
	*ff = (sin(qr)-qr*cos(qr))/(qr*qr*qr);
}
