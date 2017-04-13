/*
 * common.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef COMMON_H_
#define COMMON_H_

struct Record
{
	char * contig;             /*Contig name*/
	unsigned int position;     /*Position at contig*/
	char * referenceContext;   /*Dinucleotide base call according to reference*/
	char * callContext;        /*Dinucleotide base call according to genotype call*/
	unsigned int phredScore;   /*Probability call*/
	float methValue;           /*Methylation Value*/
	float methDev;             /*Methilation Deviation*/
	int noValue;               /*1 If there is no methylation Information*/
};

struct Bed
{
	char * contig;                   /*Contig name*/
	unsigned int start;              /*Star Position*/
	unsigned int end;                /*End Position*/
	char * extra;                    /*Extra BED Fields*/
	float meanMeth;                  /*Mean Meth*/
	float medianMeth;                /*Median Meth*/
	float stDevMeth;                 /*St Dev Methilation*/
	unsigned int cpgDinucleotides;   /*Number of CpG Dinucleotides In a given region*/
	unsigned int snps;               /*Number of snps in the bed Region*/
};



#endif /* COMMON_H_ */
