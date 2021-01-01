#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "header.h"
#include "matvec.h"

void cmalloc_(void **ptr, unsigned int *bytes)
{
	/* Allocate the memory */
	*ptr=malloc(*bytes);

	/* We will zero out the memory */
	memset(*ptr, '\0', *bytes);
	
	return;
}

void cfree_(void **ptr)
{
	/* Free the memory */
	free(*ptr);
	*ptr=NULL;
}
