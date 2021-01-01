/*
 *  output.c
 *
 *  Created by Xuehai Huang on 08/31/2010.
 *  Copyright 2010 WZU. All rights reserved.
 *
 */

/*! \file output.c
 *  \brief Output to files
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"


/**
 * \fn int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh)
 * \brief output elements and nodes data used for drawing mesh in Matlab
 * \param *elements pointer to the structure of the triangulation
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *uh pointer to the value of solution
 * \return 1: success 0: false
 */
int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh)
{
	int i;
	int Nt=elements->row;
	char *fnodes="output/nodes.dat";
	char *felements="output/elements.dat";
	char *fuh="output/uh.dat";

	FILE *outputFile;
	outputFile=fopen(fnodes, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", fnodes);
		return 0;
	}
	for(i=0;i<nodes->row;i++)
		fprintf(outputFile, "%e %e %d\n", nodes->val[i][0], nodes->val[i][1], nodes->bdFlag[i]);
	fclose(outputFile);
	
	outputFile=fopen(felements, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", felements);
		return 0;
	}
	for(i=0;i<elements->row;i++)
		fprintf(outputFile, "%d %d %d\n", elements->val[i][0]+1, elements->val[i][1]+1, elements->val[i][2]+1);
	fclose(outputFile);

	outputFile=fopen(fuh, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", fuh);
		return 0;
	}
//	for(i=0;i<nodes->row;i++)
	for (i = 0; i<uh->row; i++)
		fprintf(outputFile, "%e\n", uh->val[i]);
	fclose(outputFile);

	return 1;
}
