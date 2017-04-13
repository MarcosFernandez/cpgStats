/*
 * counts.cpp
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#include "counts.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * \brief Counts Initialization
 */
void initCounts()
{
	unsigned int i;

	for (i=0; i < 4; i=i+1)
	{
	    counts.homozygousMethylated[i] = 0;
	    counts.homozygousUndefMethylated[i] = 0;
	    counts.homozygousUnMethylated[i] = 0;

	    counts.heterozygousMethylated[i] = 0;
	    counts.heterozygousUndefMethylated[i] = 0;
	    counts.heterozygousUnMethylated[i] = 0;

	    counts.deNovoCpgMethylated[i] = 0;
	    counts.deNovoCpgUndefMethylated[i] = 0;
	    counts.deNovoCpgUnMethylated[i] = 0;

	    counts.cpgWithSnpMethylated[i] = 0;
	    counts.cpgWithSnpUndefMethylated[i] = 0;
	    counts.cpgWithSnpUnMethylated[i] = 0;
	}
}

/**
 * \brief Add a new record on a given truct
 * \param Record stats
 * \param Counter vector to Update
 */
void addCounts(struct Record * record,unsigned int * cnt)
{
	if(record->phredScore < 10)
	{
		/*Quality Under 10*/
		cnt[0] = cnt[0] + 1;
	}
	else if(record->phredScore >= 10 && record->phredScore < 20)
	{
		/*Quality between 10 and 20*/
		cnt[1] = cnt[1] + 1;
	}
	else if(record->phredScore >= 20 && record->phredScore < 30)
	{
		/*Quality between 20 and 30*/
		cnt[2] = cnt[2] + 1;
	}
	else
	{
		/*Quality below 30*/
		cnt[3] = cnt[3] + 1;
	}
}

/**
 * \brief Add record to the global stats
 * \params record -  Record Stats to add
 */
void addRecordStats(struct Record * record)
{
    /* 1. Check contig existance */
	if(record->contig == NULL)
	{
		return;
	}

	/* 2. Check Methylation Value existance */
	if(record->noValue == 1)
	{
		return;
	}


	if(strcmp(record->referenceContext, record->callContext) == 0)
	{
        /*3.1 Homozygous Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.1.1 Methylated*/
			addCounts(record,counts.homozygousUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.1.2 UnDefined Methylation*/
			addCounts(record,counts.homozygousUndefMethylated);
		}
		else
		{
			/*3.1.3 UnMethylated*/
			addCounts(record,counts.homozygousMethylated);
		}
	}
	else
	{
		/*3.2 Heterozygous Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.2.1 Methylated*/
			addCounts(record,counts.heterozygousUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.2.2 UnDefined Methylation*/
			addCounts(record,counts.heterozygousUndefMethylated);
		}
		else
		{
			/*3.2.3 UnMethylated*/
			addCounts(record,counts.heterozygousMethylated);
		}

		/*3.3 DeNovo CpGs Detectected*/
		if(strcmp(record->referenceContext, "CG") != 0 && strcmp(record->callContext, "CG") == 0 )
		{
			if (record->methValue <= 0.3)
			{
				/*3.3.1 Methylated*/
				addCounts(record,counts.deNovoCpgUnMethylated);
			}
			else if(record->methValue > 0.3 && record->methValue <= 0.7)
			{
				/*3.3.2 UnDefined Methylation*/
				addCounts(record,counts.deNovoCpgUndefMethylated);
			}
			else
			{
				/*3.3.3 UnMethylated*/
				addCounts(record,counts.deNovoCpgMethylated);
			}
		}

		/*3.4 CpG With SNPs*/
		if(strcmp(record->referenceContext, "CG") == 0 && strcmp(record->callContext, "CG") != 0 )
		{
			if (record->methValue <= 0.3)
			{
				/*3.4.1 Methylated*/
				addCounts(record,counts.cpgWithSnpUnMethylated);
			}
			else if(record->methValue > 0.3 && record->methValue <= 0.7)
			{
				/*3.4.2 UnDefined Methylation*/
				addCounts(record,counts.cpgWithSnpUndefMethylated);
			}
			else
			{
				/*3.4.3 UnMethylated*/
				addCounts(record,counts.cpgWithSnpMethylated);
			}
		}
	}
}


/**
 * \brief Get Total from a given vector
 * \param Get counter vector
 * \return Return total value
 */
unsigned int getTotalCount(unsigned int * cnt)
{
	unsigned int i = 0;
	unsigned int total = 0;

	for (i=0; i< 4; i=i+1)
	{
		total = total + cnt[i];
	}
	return total;
}



/* Get Total Homozygous Counts */
unsigned int getTotalHomozygous(){return (getTotalCount(counts.homozygousMethylated) + getTotalCount(counts.homozygousUndefMethylated) + getTotalCount(counts.homozygousUnMethylated));}

/* Get Total Homozygous Counts High Quality*/
unsigned int getTotalHomozygousHighQuality()
{
	return counts.homozygousMethylated[3] + counts.homozygousUndefMethylated[3] + counts.homozygousUnMethylated[3];
}

/* Get Total Heterozygous Counts */
unsigned int getTotalHeterozygous(){return (getTotalCount(counts.heterozygousMethylated) + getTotalCount(counts.heterozygousUndefMethylated) + getTotalCount(counts.heterozygousUnMethylated));}

/* Get Total Heterozygous Counts High Quality */
unsigned int getTotalHeterozygousHighQuality()
{
	return counts.heterozygousMethylated[3] + counts.heterozygousUndefMethylated[3] + counts.heterozygousUnMethylated[3];
}

/* Get Total Methylated */
unsigned int getTotalMethylated(){ return (getTotalCount(counts.heterozygousMethylated) + getTotalCount(counts.homozygousMethylated));}

/* Get Total Methylated High Quality*/
unsigned int getTotalMethylatedHighQuality(){ return counts.heterozygousMethylated[3] + counts.homozygousMethylated[3];}


/* Get Total Undefined Methylated */
unsigned int getTotalUndefinedMethylated(){ return (getTotalCount(counts.heterozygousUndefMethylated) + getTotalCount(counts.homozygousUndefMethylated));}

/* Get Total Undefined Methylated High Quality*/
unsigned int getTotalUndefinedMethylatedHighQuality(){ return counts.heterozygousUndefMethylated[3] + counts.homozygousUndefMethylated[3];}

/* Get Total UnMethylated */
unsigned int getTotalUnMethylated(){ return (getTotalCount(counts.homozygousUnMethylated) + getTotalCount(counts.heterozygousUnMethylated));}

/* Get Total UnMethylated High Quality*/
unsigned int getTotalUnMethylatedHighQuality(){ return counts.homozygousUnMethylated[3] + counts.heterozygousUnMethylated[3];}

/* Get Total Methylated Homozygous */
unsigned int getTotalMethylatedHomozygous(){ return getTotalCount(counts.homozygousMethylated);}

/* Get Total Methylated Homozigous High Quality */
unsigned int getTotalMethylatedHomozygousHighQuality(){ return counts.homozygousMethylated[3];}

/* Get Total Methylated Heterozygous */
unsigned int getTotalMethylatedHeterozygous(){ return getTotalCount(counts.heterozygousMethylated);}

/* Get Total Methylated Heterozygous High Quality */
unsigned int getTotalMethylatedHeterozygousHighQuality(){ return counts.heterozygousMethylated[3];}


/* Get Total Undefined Methylated Homozygous */
unsigned int getTotalUndefinedMethylatedHomozygous(){ return getTotalCount(counts.homozygousUndefMethylated);}

/* Get Total Undefined Methylated Homozigous High Quality */
unsigned int getTotalUndefinedMethylatedHomozygousHighQuality(){ return counts.homozygousUndefMethylated[3];}

/* Get Total Undefined Methylated Heterozygous */
unsigned int getTotalUndefinedMethylatedHeterozygous(){ return getTotalCount(counts.heterozygousUndefMethylated);}

/* Get Total Undefined Methylated Heterozygous High Quality */
unsigned int getTotalUndefinedMethylatedHeterozygousHighQuality(){ return counts.heterozygousUndefMethylated[3];}


/* Get Total UnMethylated Homozygous */
unsigned int getTotalUnMethylatedHomozygous(){ return getTotalCount(counts.homozygousUnMethylated);}

/* Get Total UnMethylated Homozigous High Quality */
unsigned int getTotalUnMethylatedHomozygousHighQuality(){ return counts.homozygousUnMethylated[3];}

/* Get Total UnMethylated Heterozygous */
unsigned int getTotalUnMethylatedHeterozygous(){ return getTotalCount(counts.heterozygousUnMethylated);}

/* Get Total UnMethylated Heterozygous High Quality */
unsigned int getTotalUnMethylatedHeterozygousHighQuality(){ return counts.heterozygousUnMethylated[3];}

/* Get Total number of Dinucleotides*/
unsigned int getTotalDinucleotides(){ return getTotalHomozygous() + getTotalHeterozygous();}

/*Get Total Quality Under 10 */
unsigned int getTotalQualityUnder10()
{
	return counts.homozygousMethylated[0] + counts.homozygousUndefMethylated[0] + counts.homozygousUnMethylated[0] +
	       counts.heterozygousMethylated[0] + counts.heterozygousUndefMethylated[0] + counts.heterozygousUnMethylated[0];
}

/*Get Total Quality Between 10 and 20 */
unsigned int getTotalQualityBetween10_20()
{
	return counts.homozygousMethylated[1] + counts.homozygousUndefMethylated[1] + counts.homozygousUnMethylated[1] +
		   counts.heterozygousMethylated[1] + counts.heterozygousUndefMethylated[1] + counts.heterozygousUnMethylated[1];
}

/*Get Total Quality Between 20 and 30 */
unsigned int getTotalQualityBetween20_30()
{
	return counts.homozygousMethylated[2] + counts.homozygousUndefMethylated[2] + counts.homozygousUnMethylated[2] +
		   counts.heterozygousMethylated[2] + counts.heterozygousUndefMethylated[2] + counts.heterozygousUnMethylated[2];
}

/* Get Total HighQuality */
unsigned int getTotalHighQuality()
{
	return counts.homozygousMethylated[3] + counts.homozygousUndefMethylated[3] + counts.homozygousUnMethylated[3] +
		   counts.heterozygousMethylated[3] + counts.heterozygousUndefMethylated[3] + counts.heterozygousUnMethylated[3];
}



/**
 * \brief Get Total De Novo CpGs
 */
unsigned int getTotalDeNovoCpgs()
{
	return getTotalCount(counts.deNovoCpgMethylated) + getTotalCount(counts.deNovoCpgUndefMethylated) + getTotalCount(counts.deNovoCpgUnMethylated);
}

/**
 * \brief Get Total De Novo CpG quality over 30
 */
unsigned int getTotalDeNovoCpgsHighQuality()
{
	return counts.deNovoCpgMethylated[3] + counts.deNovoCpgUndefMethylated[3] + counts.deNovoCpgUnMethylated[3];
}

/**
 * \brief Get Total de Novo
 */
unsigned int getTotalDeNovoCpgsQuality_over20()
{
	return counts.deNovoCpgMethylated[2] + counts.deNovoCpgUndefMethylated[2] + counts.deNovoCpgUnMethylated[2] +
		   counts.deNovoCpgMethylated[3] + counts.deNovoCpgUndefMethylated[3] + counts.deNovoCpgUnMethylated[3];
}


/**
 * \brief Total De Novo CpG Methylated
 */
unsigned int getTotalDeNovoCpgsMethylated(){ return getTotalCount(counts.deNovoCpgMethylated); }


/**
 * \brief Total De Novo High Quality Methylated
 */
unsigned int getTotalDeNovoCpgsMethylatedHighQuality(){	return counts.deNovoCpgMethylated[3];}

/**
 * \brief Get Total De Novo Cpgs Methylated Quality Over 20
 */
unsigned int getTotalDeNovoCpgsMethylatedQuality_over20(){	return counts.deNovoCpgMethylated[2] +   counts.deNovoCpgMethylated[3];}


/**
 * \brief Total De Novo CpG Undefined Methylated
 */
unsigned int getTotalDeNovoCpgsUndefinedMethylated(){ return getTotalCount(counts.deNovoCpgUndefMethylated);}

/**
 * \brief Total De Novo High Quality Undefined Methylated
 */
unsigned int getTotalDeNovoCpgsUndefinedMethylatedHighQuality(){ return counts.deNovoCpgUndefMethylated[3];}

/**
 * \brief Get Total De Novo Cpgs Undefined Methylated Quality Over 20
 */
unsigned int getTotalDeNovoCpgsUndefinedMethylatedQuality_over20(){ return counts.deNovoCpgUndefMethylated[2] + counts.deNovoCpgUndefMethylated[3];}


/**
 * \brief Total DeNovo CpGs UnMethylated
 */
unsigned int getTotalDeNovoCpgsUnMethylated(){ return getTotalCount(counts.deNovoCpgUnMethylated);}


/**
 * \brief Total De Novo High Quality UnMethylated
 */
unsigned int getTotalDeNovoCpgsUnMethylatedHighQuality(){ return counts.deNovoCpgUnMethylated[3]; }

/**
 * \brief Get Total De Novo Cpgs UnMethylated Quality Over 20
 */
unsigned int getTotalDeNovoCpgsUnMethylatedQuality_over20(){ return counts.deNovoCpgUnMethylated[2] + counts.deNovoCpgUnMethylated[3];}


/**
 * \brief Total CpG with SNP Called
 */
unsigned int getTotalCpgSnp()
{
	return getTotalCount(counts.cpgWithSnpMethylated) + getTotalCount(counts.cpgWithSnpUndefMethylated) + getTotalCount(counts.cpgWithSnpUnMethylated);
}

/**
 * \brief Total CpG with SNP Called High Quality
 */
unsigned int getTotalCpgSnpHighQuality()
{
	return  counts.deNovoCpgMethylated[3] + counts.deNovoCpgUndefMethylated[3] + counts.deNovoCpgUnMethylated[3];
}

/**
 * \brief Total CpG with SNP Called Quality Over 20
 */
unsigned int getTotalCpgSnpQuality_over20()
{
	return  counts.deNovoCpgMethylated[2] + counts.deNovoCpgUndefMethylated[2] + counts.deNovoCpgUnMethylated[2] +
			counts.deNovoCpgMethylated[3] + counts.deNovoCpgUndefMethylated[3] + counts.deNovoCpgUnMethylated[3];
}


/**
 * \brief Total CpGs With SNPs Methylated
 */
unsigned int getTotalCpgSnpMethylated(){ return getTotalCount(counts.cpgWithSnpMethylated);}


/**
 * \brief Total CpG with SNP Called Methylated High Quality
 */
unsigned int getTotalCpgSnpMethylatedHighQuality(){ return counts.cpgWithSnpMethylated[3]; }


/**
 * \brief Total CpG with SNP Called Methylated Quality over 20
 */
unsigned int getTotalCpgSnpMethylatedQuality_over20(){ return counts.cpgWithSnpMethylated[2] + counts.cpgWithSnpMethylated[3]; }




/**
 * \brief Total CpGs With SNPs Undefined Methylated
 */
unsigned int getTotalCpgSnpUndefinedMethylated(){ return getTotalCount(counts.cpgWithSnpUndefMethylated);}


/**
 * \brief Total CpGs With SNPs Undefined Methylated High Quality
 */
unsigned int getTotalCpgSnpUndefinedMethylatedHighQuality(){ return counts.cpgWithSnpUndefMethylated[3]; }


/**
 * \brief  Total CpGs With SNPs Undefined Methylated Quality Over 20
 */
unsigned int getTotalCpgSnpUndefinedMethylatedQuality_over20(){	return counts.cpgWithSnpUndefMethylated[2] + counts.cpgWithSnpUndefMethylated[3];}



/**
 * \brief Get Total Cpg Snp UnMethylated
 */
unsigned int getTotalCpgSnpUnMethylated(){	return getTotalCount(counts.cpgWithSnpUnMethylated);}

/**
 * \brief Get Total Cpg Snp UnMethylated High Quality
 */
unsigned int getTotalCpgSnpUnMethylatedHighQuality(){ return counts.cpgWithSnpUnMethylated[3];}

/**
 * \brief Get Total Cpg Snp UnMethylated Quality Over 20
 */
unsigned int getTotalCpgSnpUnMethylatedQuality_over20(){ return counts.cpgWithSnpUnMethylated[2] + counts.cpgWithSnpUnMethylated[3];}






/**
 * \brief get percentage value
 * \param concept to calculate the percentage
 * \param total Value from which estimate the percentage
 */
float getPercentage(unsigned int concept, unsigned int totalValue)
{
	return ((float)concept/totalValue)*100;
}

/**
 * \brief Print Counters
 */
void printCounts()
{
	printf("CpG Stats \n");
	printf("Total Dinucleotides: \t %i \t (%.2f %%) \n", getTotalDinucleotides(),getPercentage(getTotalDinucleotides(), getTotalDinucleotides()));
	printf("\t Total High Quality: \t %i \t (%.2f %%) \n", getTotalHighQuality(),getPercentage(getTotalHighQuality(), getTotalDinucleotides()));
	printf("\t Total Quality (20-30): \t %i \t (%.2f %%) \n", getTotalQualityBetween20_30(),getPercentage(getTotalQualityBetween20_30(), getTotalDinucleotides()));
	printf("\t Total Quality (10-20): \t %i \t (%.2f %%) \n", getTotalQualityBetween10_20(),getPercentage(getTotalQualityBetween10_20(), getTotalDinucleotides()));
	printf("\t Total Quality Under 10: \t %i \t (%.2f %%) \n", getTotalQualityUnder10(),getPercentage(getTotalQualityUnder10(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Homozygous: \t %i \t (%.2f %%) \n",getTotalHomozygous(),getPercentage(getTotalHomozygous(), getTotalDinucleotides()));
	printf("\t Total Homozygous High Quality: \t %i \t (%.2f %%) \n",getTotalHomozygousHighQuality(),getPercentage(getTotalHomozygousHighQuality(), getTotalDinucleotides()));
	printf("Total Heterozygous: \t %i \t (%.2f %%) \n",getTotalHeterozygous(),getPercentage(getTotalHeterozygous(), getTotalDinucleotides()));
	printf("\t Total Heterozygous High Quality: \t %i \t (%.2f %%) \n",getTotalHeterozygousHighQuality(),getPercentage(getTotalHeterozygousHighQuality(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Methylated: \t %i \t (%.2f %%) \n", getTotalMethylated(),getPercentage(getTotalMethylated(), getTotalDinucleotides()));
	printf("\t Total Methylated High Quality: \t %i \t (%.2f %%) \n", getTotalMethylatedHighQuality(),getPercentage(getTotalMethylatedHighQuality(), getTotalDinucleotides()));
    printf("Total Undefined Methylated: \t %i \t (%.2f %%) \n", getTotalUndefinedMethylated(),getPercentage(getTotalUndefinedMethylated(), getTotalDinucleotides()));
	printf("\t Total Undefined Methylated High Quality: \t %i \t (%.2f %%) \n", getTotalUndefinedMethylatedHighQuality(),getPercentage(getTotalUndefinedMethylatedHighQuality(), getTotalDinucleotides()));
	printf("Total UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylated(),getPercentage(getTotalUnMethylated(), getTotalDinucleotides()));
	printf("\t Total Unmethylated High Quality: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHighQuality(),getPercentage(getTotalUnMethylatedHighQuality(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Homozygous Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedHomozygous(),getPercentage(getTotalMethylatedHomozygous(),getTotalMethylated()));
	printf("\t Total Homozygous Methylated High Quality: \t %i \t (%.2f %%)  \n", getTotalMethylatedHomozygousHighQuality(),getPercentage(getTotalMethylatedHomozygousHighQuality(),getTotalMethylated()));
	printf("Total Heterozygous Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedHeterozygous(),getPercentage(getTotalMethylatedHeterozygous(),getTotalMethylated()));
	printf("\t Total Heterozygous Methylated High Quality: \t %i \t (%.2f %%)  \n", getTotalMethylatedHeterozygousHighQuality(),getPercentage(getTotalMethylatedHeterozygousHighQuality(),getTotalMethylated()));

	printf("\n");
	printf("Total Homozygous Undefined Methylated: \t %i \t (%.2f %%)  \n", getTotalUndefinedMethylatedHomozygous(),getPercentage(getTotalUndefinedMethylatedHomozygous(),getTotalUndefinedMethylated()));
	printf("\t Total Homozygous Undefined Methylated High Quality: \t %i \t (%.2f %%)  \n", getTotalUndefinedMethylatedHomozygousHighQuality(),getPercentage(getTotalUndefinedMethylatedHomozygousHighQuality(),getTotalUndefinedMethylated()));
	printf("Total Heterozygous Undefined Methylated: \t %i \t (%.2f %%)  \n", getTotalUndefinedMethylatedHeterozygous(),getPercentage(getTotalUndefinedMethylatedHeterozygous(),getTotalUndefinedMethylated()));
	printf("\t Total Heterozygous Undefined Methylated High Quality: \t %i \t (%.2f %%)  \n", getTotalUndefinedMethylatedHeterozygousHighQuality(),getPercentage(getTotalUndefinedMethylatedHeterozygousHighQuality(),getTotalUndefinedMethylated()));

	printf("\n");
	printf("Total Homozygous UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHomozygous(),getPercentage(getTotalUnMethylatedHomozygous(),getTotalUnMethylated()));
	printf("\t Total Homozygous UnMethylated High Quality: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHomozygousHighQuality(),getPercentage(getTotalUnMethylatedHomozygousHighQuality(),getTotalUnMethylated()));
	printf("Total Heterozygous UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHeterozygous(),getPercentage(getTotalUnMethylatedHeterozygous(),getTotalUnMethylated()));
	printf("\t Total Heterozygous UnMethylated High Quality: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHeterozygousHighQuality(),getPercentage(getTotalUnMethylatedHeterozygousHighQuality(),getTotalUnMethylated()));

	printf("\n");
	printf("Total DeNovo CpGs: \t %i \t (%.2f %%) \n", getTotalDeNovoCpgs(),getPercentage(getTotalDeNovoCpgs(),getTotalDinucleotides()));
	printf("\t Total DeNovo CpGs Quality > 20: \t %i \t (%.2f %%) \n", getTotalDeNovoCpgsQuality_over20(),getPercentage(getTotalDeNovoCpgsQuality_over20(),getTotalDinucleotides()));
	printf("\t Total DeNovo CpGs High Quality: \t %i \t (%.2f %%) \n", getTotalDeNovoCpgsHighQuality(),getPercentage(getTotalDeNovoCpgsMethylatedHighQuality(),getTotalDinucleotides()));

	printf("\t Total DeNovo CpGs Methylated: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsMethylated(),getPercentage(getTotalDeNovoCpgsMethylated(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsMethylatedQuality_over20(),getPercentage(getTotalDeNovoCpgsMethylatedQuality_over20(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs Methylated High Quality: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsMethylatedHighQuality(),getPercentage(getTotalDeNovoCpgsMethylatedHighQuality(),getTotalDeNovoCpgs()));

	printf("\t Total DeNovo CpGs Undefined Methylated: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUndefinedMethylated(),getPercentage(getTotalDeNovoCpgsUndefinedMethylated(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs Undefined Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUndefinedMethylatedQuality_over20(),getPercentage(getTotalDeNovoCpgsUndefinedMethylatedQuality_over20(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs Undefined Methylated High Quality: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUndefinedMethylatedHighQuality(),getPercentage(getTotalDeNovoCpgsUndefinedMethylatedHighQuality(),getTotalDeNovoCpgs()));

	printf("\t Total DeNovo CpGs UnMethylated: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUnMethylated(),getPercentage(getTotalDeNovoCpgsUnMethylated(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs UnMethylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUnMethylatedQuality_over20(),getPercentage(getTotalDeNovoCpgsUnMethylatedQuality_over20(),getTotalDeNovoCpgs()));
	printf("\t\t Total DeNovo CpGs UnMethylated High Quality: \t %i \t (%.2f %%) \n",getTotalDeNovoCpgsUnMethylatedHighQuality(),getPercentage(getTotalDeNovoCpgsUnMethylatedHighQuality(),getTotalDeNovoCpgs()));

	printf("\n");
	printf("Total CpG with SNP Called: \t %i \t (%.2f %%) \n", getTotalCpgSnp(),getPercentage(getTotalCpgSnp(),getTotalDinucleotides()));
	printf("\t Total CpG with SNP Called Quality > 20: \t %i \t (%.2f %%) \n", getTotalCpgSnpQuality_over20(),getPercentage(getTotalCpgSnpQuality_over20(),getTotalDinucleotides()));
	printf("\t Total CpG with SNP Called High Quality: \t %i \t (%.2f %%) \n", getTotalCpgSnpHighQuality(),getPercentage(getTotalCpgSnpMethylatedHighQuality(),getTotalDinucleotides()));

	printf("\t Total CpG with SNP Called Methylated: \t %i \t (%.2f %%) \n",getTotalCpgSnpMethylated(),getPercentage(getTotalCpgSnpMethylated(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalCpgSnpMethylatedQuality_over20(),getPercentage(getTotalCpgSnpMethylatedQuality_over20(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called Methylated High Quality: \t %i \t (%.2f %%) \n",getTotalCpgSnpMethylatedHighQuality(),getPercentage(getTotalCpgSnpMethylatedHighQuality(),getTotalCpgSnp()));

	printf("\t Total CpG with SNP Called Undefined Methylated: \t %i \t (%.2f %%) \n",getTotalCpgSnpUndefinedMethylated(),getPercentage(getTotalCpgSnpUndefinedMethylated(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called Undefined Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalCpgSnpUndefinedMethylatedQuality_over20(),getPercentage(getTotalCpgSnpUndefinedMethylatedQuality_over20(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called Undefined Methylated High Quality: \t %i \t (%.2f %%) \n",getTotalCpgSnpUndefinedMethylatedHighQuality(),getPercentage(getTotalCpgSnpUndefinedMethylatedHighQuality(),getTotalCpgSnp()));

	printf("\t Total CpG with SNP Called UnMethylated: \t %i \t (%.2f %%) \n",getTotalCpgSnpUnMethylated(),getPercentage(getTotalCpgSnpUnMethylated(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called UnMethylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalCpgSnpUnMethylatedQuality_over20(),getPercentage(getTotalCpgSnpUnMethylatedQuality_over20(),getTotalCpgSnp()));
	printf("\t\t Total CpG with SNP Called UnMethylated High Quality: \t %i \t (%.2f %%) \n",getTotalCpgSnpUnMethylatedHighQuality(),getPercentage(getTotalCpgSnpUnMethylatedHighQuality(),getTotalCpgSnp()));

}

/**
 * \brief Print to JSON file
 * \param json file to store the data
 */
void saveCounts(char * fileName)
{
	FILE *fp;

	fp = fopen(fileName, "w");

	fprintf(fp,"{\n");

	fprintf(fp,"  \"TotalDinucleotides\":%i,\n", getTotalDinucleotides());
	fprintf(fp,"  \"TotalHighQuality\":%i,\n", getTotalHighQuality());
	fprintf(fp,"  \"TotalQuality_20_30\":%i,\n", getTotalQualityBetween20_30());
	fprintf(fp,"  \"TotalQuality_10_20\":%i,\n", getTotalQualityBetween10_20());
	fprintf(fp,"  \"TotalQualityUnder_10\":%i,\n", getTotalQualityUnder10());

	fprintf(fp,"  \"TotalHomozygous\":%i,\n",getTotalHomozygous());
	fprintf(fp,"  \"TotalHomozygousHighQuality\":%i,\n",getTotalHomozygousHighQuality());
	fprintf(fp,"  \"TotalHeterozygous\":%i,\n",getTotalHeterozygous());
	fprintf(fp,"  \"TotalHeterozygousHighQuality\":%i,\n",getTotalHeterozygousHighQuality());

	fprintf(fp,"  \"TotalMethylated\":%i,\n", getTotalMethylated());
	fprintf(fp,"  \"TotalMethylatedHighQuality\":%i,\n", getTotalMethylatedHighQuality());
	fprintf(fp,"  \"TotalUndefinedMethylated\":%i,\n", getTotalUndefinedMethylated());
	fprintf(fp,"  \"TotalUndefinedMethylatedHighQuality\":%i,\n", getTotalUndefinedMethylatedHighQuality());
	fprintf(fp,"  \"TotalUnMethylated\":%i,\n", getTotalUnMethylated());
	fprintf(fp,"  \"TotalUnmethylatedHighQuality\":%i,\n", getTotalUnMethylatedHighQuality());

	fprintf(fp,"  \"TotalHomozygousMethylated\":%i,\n", getTotalMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousMethylatedHighQuality\":%i,\n", getTotalMethylatedHomozygousHighQuality());
	fprintf(fp,"  \"TotalHeterozygousMethylated\":%i,\n", getTotalMethylatedHeterozygous());
	fprintf(fp,"  \"TotalHeterozygousMethylatedHighQuality\":%i,\n", getTotalMethylatedHeterozygousHighQuality());

	fprintf(fp,"  \"TotalHomozygousUndefinedMethylated\":%i,\n", getTotalUndefinedMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousUndefinedMethylatedHighQuality\":%i,\n", getTotalUndefinedMethylatedHomozygousHighQuality());
	fprintf(fp,"  \"TotalHeterozygousUndefinedMethylated\":%i,\n", getTotalUndefinedMethylatedHeterozygous());
	fprintf(fp,"  \"TotalHeterozygousUndefinedMethylatedHighQuality\":%i,\n", getTotalUndefinedMethylatedHeterozygousHighQuality());

	fprintf(fp,"  \"TotalHomozygousUnMethylated\":%i,\n", getTotalUnMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousUnMethylatedHighQuality\":%i,\n", getTotalUnMethylatedHomozygousHighQuality());
	fprintf(fp,"  \"TotalHeterozygousUnMethylated\":%i,\n", getTotalUnMethylatedHeterozygous());
	fprintf(fp,"  \"TotalHeterozygousUnMethylatedHighQuality\":%i,\n", getTotalUnMethylatedHeterozygousHighQuality());

	fprintf(fp,"  \"TotalDeNovoCpGs\":%i,\n", getTotalDeNovoCpgs());
	fprintf(fp,"  \"TotalDeNovoCpGsQualityO20\":%i,\n", getTotalDeNovoCpgsQuality_over20());
	fprintf(fp,"  \"TotalDeNovoCpGsHighQuality\":%i,\n", getTotalDeNovoCpgsHighQuality());

	fprintf(fp,"  \"TotalDeNovoCpGsMethylated\":%i,\n",getTotalDeNovoCpgsMethylated());
	fprintf(fp,"  \"TotalDeNovoCpGsMethylatedQualityO20\":%i,\n",getTotalDeNovoCpgsMethylatedQuality_over20());
	fprintf(fp,"  \"TotalDeNovoCpGsMethylatedHighQuality\":%i,\n",getTotalDeNovoCpgsMethylatedHighQuality());

	fprintf(fp,"  \"TotalDeNovoCpGsUndefinedMethylated\":%i,\n",getTotalDeNovoCpgsUndefinedMethylated());
	fprintf(fp,"  \"TotalDeNovoCpGsUndefinedMethylatedQualityO20\":%i,\n",getTotalDeNovoCpgsUndefinedMethylatedQuality_over20());
	fprintf(fp,"  \"TotalDeNovoCpGsUndefinedMethylatedHighQuality\":%i,\n",getTotalDeNovoCpgsUndefinedMethylatedHighQuality());

	fprintf(fp,"  \"TotalDeNovoCpGsUnMethylated\":%i,\n",getTotalDeNovoCpgsUnMethylated());
	fprintf(fp,"  \"TotalDeNovoCpGsUnMethylatedQualityO20\":%i,\n",getTotalDeNovoCpgsUnMethylatedQuality_over20());
	fprintf(fp,"  \"TotalDeNovoCpGsUnMethylatedHighQuality\":%i,\n",getTotalDeNovoCpgsUnMethylatedHighQuality());

	fprintf(fp,"  \"TotalCpGwithSNPCalled\":%i,\n", getTotalCpgSnp());
	fprintf(fp,"  \"TotalCpGwithSNPCalledQualityO20\":%i,\n", getTotalCpgSnpQuality_over20());
	fprintf(fp,"  \"TotalCpGwithSNPCalledHighQuality\":%i,\n", getTotalCpgSnpHighQuality());

	fprintf(fp,"  \"TotalCpGwithSNPCalledMethylated\":%i,\n",getTotalCpgSnpMethylated());
	fprintf(fp,"  \"TotalCpGwithSNPCalledMethylatedQualityO20\":%i,\n",getTotalCpgSnpMethylatedQuality_over20());
	fprintf(fp,"  \"TotalCpGwithSNPCalledMethylatedHighQuality\":%i,\n",getTotalCpgSnpMethylatedHighQuality());

	fprintf(fp,"  \"TotalCpGwithSNPCalledUndefinedMethylated\":%i,\n",getTotalCpgSnpUndefinedMethylated());
	fprintf(fp,"  \"TotalCpGwithSNPCalledUndefinedMethylatedQualityO20:\":%i,\n",getTotalCpgSnpUndefinedMethylatedQuality_over20());
	fprintf(fp,"  \"TotalCpGwithSNPCalledUndefinedMethylatedHighQuality\":%i,\n",getTotalCpgSnpUndefinedMethylatedHighQuality());

	fprintf(fp,"  \"TotalCpGwithSNPCalledUnMethylated\":%i,\n",getTotalCpgSnpUnMethylated());
	fprintf(fp,"  \"TotalCpGwithSNPCalledUnMethylatedQualityO20\":%i,\n",getTotalCpgSnpUnMethylatedQuality_over20());
	fprintf(fp,"  \"TotalCpGwithSNPCalledUnMethylatedHighQuality\":%i\n",getTotalCpgSnpUnMethylatedHighQuality());

	fprintf(fp,"}\n");

	fclose(fp);
}










