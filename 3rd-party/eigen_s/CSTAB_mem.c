#include <stdio.h>


void
get_delta_(void *a1, void *a2, int *i)
{
	long	b1 = (long)a1;
	long	b2 = (long)a2;
	long	delta;

	delta = b1-b2;
	*i = (int) delta;
}


