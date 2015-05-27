/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate 
# students at Oregon State University; and Dr. Nik Grünwald, an employee of 
# USDA-ARS.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee, 
# and without a written agreement is hereby granted, provided that the statement
# above is incorporated into the material, giving appropriate attribution to the
# authors.
#
# Permission to incorporate this software into commercial products may be
# obtained by contacting USDA ARS and OREGON STATE UNIVERSITY Office for 
# Commercialization and Corporate Development.
#
# The software program and documentation are supplied "as is", without any
# accompanying services from the USDA or the University. USDA ARS or the 
# University do not warrant that the operation of the program will be 
# uninterrupted or error-free. The end-user understands that the program was 
# developed for research purposes and is advised not to rely exclusively on the 
# program for any reason.
#
# IN NO EVENT SHALL USDA ARS OR OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY 
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
# LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
# EVEN IF THE OREGON STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF 
# SUCH DAMAGE. USDA ARS OR OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY 
# WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
# BASIS, AND USDA ARS AND OREGON STATE UNIVERSITY HAVE NO OBLIGATIONS TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
#
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <Rinternals.h>
#include <R.h>
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A slightly faster method of permuting alleles at a locus. 

Inputs:
	locus - The matrix to be permuted. Used for reference.
	alleles - a vector of integers from 0 to n alleles indicating the matrix cols
	# allele_freq - 1/ploidy
	ploidy - self-explainitory

Outputs;
	Rout - the permuted matrix. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP permute_shuff(SEXP locus, SEXP alleles, SEXP ploidy)
{
	int rows;
	int cols;
	int i;
	int j;
	int count = 0;
	int* ploid;
	int p;
	int miss = 0;
	SEXP Rout;
	SEXP Rdim;
	Rdim = getAttrib(locus, R_DimSymbol);
	rows = INTEGER(Rdim)[0]; 
	cols = INTEGER(Rdim)[1]; 
	PROTECT(Rout = allocMatrix(INTSXP, rows, cols));
	alleles = coerceVector(alleles, INTSXP);
	ploidy = coerceVector(ploidy, INTSXP);
	ploid = INTEGER(ploidy);
	int* inmat = INTEGER(locus);
	int* outmat = INTEGER(Rout);
	int* alle = INTEGER(alleles);
	for(i = 0; i < rows; i++)
	{
		// loop through all columns first and initialize
		for(j = 0; j < cols; j++) 
		{
			if (inmat[i + j*rows] == NA_INTEGER) // skip missing values
			{
				outmat[i + j*rows] = NA_INTEGER; 
				miss = 1;
			}
			else // initialize to zero
			{
				outmat[i + j*rows] = 0;
			}
		}
		if (miss == 1)
		{
			miss = 0;
		}
		else
		{
			// loop over the observed ploidy and assign that many alleles.
			for(p = 0; p < ploid[i]; p++)
			{
				outmat[i + alle[count++]*rows] += 1;
			}
		}
	}
	UNPROTECT(1);
	return Rout;
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Method for expanding indices for bootstrapping. Only slightly faster than R
version, but seems to scale better.

Inputs:
	indices - cumulative sum of the number of alleles at each locus.
	length - number of loci

Outputs;
	res - a list the same length as the number of loci with continuous numbers. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP expand_indices(SEXP indices, SEXP length) {
	SEXP res;
	SEXP tempvec;
	int rows;
	int i;
	int j;
	int *ind = INTEGER(indices);
	int max;
	int min = 1;
	rows = INTEGER(length)[0];
	PROTECT(res = allocVector(VECSXP, rows));
	for (i = 0; i < rows; i++)
	{
		max = ind[i];
		int veclength = 1 + max - min;
		PROTECT(tempvec = allocVector(INTSXP, veclength));
		for (j = 0; j < veclength ; j++)
		{
			INTEGER(tempvec)[j] = min + j;
		}
		SET_VECTOR_ELT(res, i, tempvec);
		UNPROTECT(1);
		min = ind[i] + 1;
	}
	UNPROTECT(1);
	return res;
}
