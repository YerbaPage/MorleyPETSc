/*
 *  sort.c
 *
 *  Created by Xuehai Huang on 08/26/2010.
 *  Copyright 2010 WZU. All rights reserved.
 *
 */

/*! \file sort.c
 *  \brief sort algorithm
 */

#include <stdlib.h>
#include <math.h>
#include "header.h"


/**
 * \fn int partitions(double *a,int low,int high)
 * \brief partitions used in quicksort algorithm
          Reorder the list so that all elements with values less than the pivot come before the pivot, 
		  while all elements with values greater than the pivot come after it (equal values can go either way). 
		  After this partitioning, the pivot is in its final position. This is called the partition operation
 * \param *a pointer to the double array
 * \param low the starting index
 * \param high the ending index   
 * \return pivot index 
 */
int partitions(double *a,int low,int high)
{
	double pivotkey=a[low];
	while(low<high)
	{
		while(low<high && a[high]>=pivotkey)
			--high;
		a[low]=a[high];
		while(low<high && a[low]<=pivotkey)
			++low;
		a[high]=a[low];
	}
	a[low]=pivotkey;
	
	return low;
}

/**
 * \fn int partitionsIndex(double *a,int *index,int low,int high)
 * \brief partitions used in quicksort algorithm with index
          Reorder the list so that all elements with values less than the pivot come before the pivot, 
		  while all elements with values greater than the pivot come after it (equal values can go either way). 
		  After this partitioning, the pivot is in its final position. This is called the partition operation
 * \param *a pointer to the double array
 * \param *index pointer to the index of 'a'
 * \param low the starting index
 * \param high the ending index   
 * \return pivot index 
 */
int partitionsIndex(double *a,int *index,int low,int high)
{
	double pivotkey=a[low];
	int indexkey=index[low];
	while(low<high)
	{
		while(low<high && a[high]>=pivotkey)
			--high;
		a[low]=a[high];
		index[low]=index[high];
		while(low<high && a[low]<=pivotkey)
			++low;
		a[high]=a[low];
		index[high]=index[low];
	}
	a[low]=pivotkey;
	index[low]=indexkey;
	
	return low;
}

/**
 * \fn void quicksort(double *a,int low,int high)
 * \brief quicksort algorithm, sort a in ascending order
 * \param *a pointer to the double array
 * \param low the starting index
 * \param high the ending index   
 * \return void 
 */
void quicksort(double *a,int low,int high)
{
	int pivottag;
	if(low<high)
	{ 
		// recursion
		pivottag=partitions(a,low,high);
		quicksort(a,low,pivottag-1);
		quicksort(a,pivottag+1,high);
	}
}

/**
 * \fn void quicksortIndex(double *a,int *index,int low,int high)
 * \brief quicksort algorithm
          reorder the index of 'a'(double type) so that 'a' is ascending in such order  
          'low' and 'high' are usually set to be 0 and n-1,respectively,where n is the 
          length of 'a'. 'index' should be initialized in the nature order and it has the
          same length as 'a'.
 * \param *a pointer to the double array
 * \param *index pointer to the index of 'a'
 * \param low the starting index
 * \param high the ending index   
 * \return void 
 */
void quicksortIndex(double *a,int *index,int low,int high)
{
	int pivottag;
	if(low<high)
	{ 
		// recursion
		pivottag=partitionsIndex(a,index,low,high);
		quicksortIndex(a,index,low,pivottag-1);
		quicksortIndex(a,index,pivottag+1,high);
	}
}