//util,cpp Mike Williams 2-05
#include"util.h"

double abs(double x)
{
	if(x<0)
		return -1*x;
	else
		return x;
}

double max(double a, double b)
{
	if(a>b)
		return a;
	return b;
}

double min(double a, double b)
{
	if(a<b)
		return a;
	return b;
}
