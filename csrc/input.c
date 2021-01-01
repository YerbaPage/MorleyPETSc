/*
 *		input.c 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 04/03/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file input.c
 *  \brief Read input parameters
 */

#include <stdio.h>
#include <string.h>

#include "header.h"
#include "matvec.h"

/**
 * \fn int read_Input_data(char *filenm, Input_data *Input)
 * \brief Read input parameters from disk file
 * \param *filenm char filename for input file
 * \param *Input input data structure
 * \return 1 if succeed, 0 if fail
 */
int read_Input_data(char *filenm, Input_data *Input)
{
	int val;
	int ibuff; double dbuff;
	char buffer[500];
	FILE *fp;
	int ret = 1;
	
	// set initial/default values in Input
	Input->print_level = 0;
	Input->problem_num = 1;
	Input->domain_num = 1;
	Input->itsolver_type = 1;
	Input->itsolver_tol = 1e-8;
	Input->itsolver_maxit = 500;
	Input->restart = 10;
	
	Input->precond_type = 1;
	Input->precond_scale[0] = 1.0;
	Input->precond_scale[1] = 1.0;
	Input->smoother = 1;
	Input->schwarz_type = 1;
	Input->smooth_iter = 1;

	Input->AMG_levels = 20;
	Input->MG_tol = 1e-8;
	Input->MG_maxit = 3;
	Input->MG_smoother = GS;
	Input->MG_smooth_iter = 1;
	Input->AMG_coarsening_type = 1;	
	Input->AMG_interpolation_type = 1;
	Input->AMG_coarse_dof = 50;	

	Input->AMG_strong_threshold = 0.5; 
	Input->AMG_truncation_threshold = 0.4;
	Input->AMG_max_row_sum = 0.6;
	
	Input->AFEM_tol = 0.001;
	Input->AFEM_maxit = 10;
	Input->AFEM_mark_threshold = 0.5;
	
	Input->nu = 0.3;
	Input->lambda = 0.3;
	Input->mu = 0.35;
	Input->t = 0;
	Input->glevelNum = 6;
	Input->CglevelNum = 4;
	Input->FglevelNum = 6;
	Input->cDOP = 2;

	Input->variationalform_type = 1;
	Input->stress_fem_type = 1;
	
	// use default input parameters
	if(filenm==NULL) return(0);
	
	fp = fopen(filenm,"r");
	if(fp==NULL){
		printf("Could not open file...\n");
		return(0);
	}
	
	while(ret)
	{
		val = fscanf(fp,"%s",buffer);
		if(val==EOF) break;
		if(val!=1){ ret=0; break; }
		if(buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
			fgets(buffer,500,fp); // skip rest of line
			continue;
		}
		
		// match keyword and scan for value
		if (strcmp(buffer,"problem_num")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->problem_num=ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer, "domain_num") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->domain_num = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer,"print_level")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->print_level = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if(strcmp(buffer,"itsolver_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 1; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 1; break; }
			Input->itsolver_type = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"itsolver_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 1; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 1; break; }
			Input->itsolver_tol = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if(strcmp(buffer,"itsolver_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 1; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 1; break; }
			Input->itsolver_maxit = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if(strcmp(buffer,"itsolver_restart")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 1; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 1; break; }
			Input->restart = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer, "precond_type") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->precond_type = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "precond_scale1") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%lf", &dbuff);
			if (val != 1) { ret = 0; break; }
			Input->precond_scale[0] = dbuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "precond_scale2") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%lf", &dbuff);
			if (val != 1) { ret = 0; break; }
			Input->precond_scale[1] = dbuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "smoother") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%s", buffer);
			if (val != 1) { ret = 0; break; }

			if (strcmp(buffer, "JACOBI") == 0)
				Input->smoother = JACOBI;
			else if (strcmp(buffer, "GS") == 0)
				Input->smoother = GS;
			else if (strcmp(buffer, "SGS") == 0)
				Input->smoother = SGS;
			else if (strcmp(buffer, "MSWZ") == 0)
				Input->smoother = MSWZ;
			else if (strcmp(buffer, "ASWZ") == 0)
				Input->smoother = ASWZ;
			else if (strcmp(buffer, "SMSWZ") == 0)
				Input->smoother = SMSWZ;
			else
			{
				ret = 0; break;
			}
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "schwarz_type") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->schwarz_type = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "smooth_iter") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->smooth_iter = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer,"MG_tol")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->MG_tol = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if(strcmp(buffer,"MG_maxit")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->MG_maxit = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer, "MG_smoother") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%s", buffer);
			if (val != 1) { ret = 0; break; }

			if (strcmp(buffer, "JACOBI") == 0)
				Input->MG_smoother = JACOBI;
			else if (strcmp(buffer, "GS") == 0)
				Input->MG_smoother = GS;
			else if (strcmp(buffer, "SGS") == 0)
				Input->MG_smoother = SGS;
			else
			{
				ret = 0; break;
			}
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "MG_smooth_iter") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->MG_smooth_iter = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "AMG_levels") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->AMG_levels = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if(strcmp(buffer,"AMG_coarse_dof")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_coarse_dof = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if(strcmp(buffer,"AMG_coarsening_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_coarsening_type = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if(strcmp(buffer,"AMG_interpolation_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_interpolation_type = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
				
		else if (strcmp(buffer,"AMG_strong_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_strong_threshold = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_truncation_threshold")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_truncation_threshold = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"AMG_max_row_sum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->AMG_max_row_sum = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer, "AFEM_tol") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%lf", &dbuff);
			if (val != 1) { ret = 0; break; }
			Input->AFEM_tol = dbuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "AFEM_maxit") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->AFEM_maxit = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "AFEM_mark_threshold") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%lf", &dbuff);
			if (val != 1) { ret = 0; break; }
			Input->AFEM_mark_threshold = dbuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer,"glevelNum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->glevelNum = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"CglevelNum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->CglevelNum = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"FglevelNum")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->FglevelNum = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"cDOP")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->cDOP = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"rDOP")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->rDOP = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"dop1")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->dop1 = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"dop2")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->dop2 = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"dop3")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->dop3 = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"dop4")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if(val!=1) { ret = 0; break; }
			Input->dop4 = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer, "variationalform_type") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->variationalform_type = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer, "stress_fem_type") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%d", &ibuff);
			if (val != 1) { ret = 0; break; }
			Input->stress_fem_type = ibuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer,"nu")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->nu = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"lambda")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->lambda = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"mu")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->mu = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer, "t") == 0)
		{
			val = fscanf(fp, "%s", buffer);
			if (val != 1 || strcmp(buffer, "=") != 0) {
				ret = 0; break;
			}
			val = fscanf(fp, "%lf", &dbuff);
			if (val != 1) { ret = 0; break; }
			Input->t = dbuff;
			fgets(buffer, 500, fp); // skip rest of line
		}

		else if (strcmp(buffer,"alpha1")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->alpha1 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"alpha2")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->alpha2 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"alpha3")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->alpha3 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"beta1")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->beta1 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"beta2")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->beta2 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}

		else if (strcmp(buffer,"beta3")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if(val!=1 || strcmp(buffer,"=")!=0) {
				ret = 0; break;
			}
			val = fscanf(fp,"%lf",&dbuff);
			if(val!=1) { ret = 0; break; }
			Input->beta3 = dbuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else
		{
			printf("Unknown Input keyword...\n");
			ret = 0;
			break;
		}		
	}
	
	if(ret!=1) {
		printf("Read_input_data: an error occured while reading the file.\n");
		return(ret);
	}
	
	// some sanity checks
	if(Input->problem_num<0 || Input->print_level<0 ||
		 Input->itsolver_tol<=0 || Input->itsolver_maxit<=0 ||
		 Input->AMG_levels<0 || Input->MG_tol<0 || Input->MG_maxit<0 ||
		 Input->AMG_coarsening_type<0 || Input->AMG_interpolation_type>100 ||
		 Input->MG_smoother<0 || Input->AMG_strong_threshold<0 || 
		 Input->AMG_truncation_threshold<0 || Input->AMG_max_row_sum<0||		 
		 Input->MG_smooth_iter<0 ) ret = 0;
	
	return(ret);
	
}
