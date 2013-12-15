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
	allele_freq - 1/ploidy
	ploidy - self-explainitory

Outputs;
	Rout - the permuted matrix. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
SEXP permute_shuff(SEXP locus, SEXP alleles, SEXP allele_freq, SEXP ploidy)
{
	int rows, cols, i, j, count = 0, ploid, p, miss = 0;
	SEXP Rout;
	SEXP Rdim;
	Rdim = getAttrib(locus, R_DimSymbol);
	rows = INTEGER(Rdim)[0]; 
	cols = INTEGER(Rdim)[1]; 
	PROTECT(Rout = allocMatrix(REALSXP, rows, cols));
	alleles = coerceVector(alleles, INTSXP);
	ploidy = coerceVector(ploidy, INTSXP);
	ploid = INTEGER(ploidy)[0];
	double *loc = REAL(locus), *outmat = REAL(Rout), *afreq = REAL(allele_freq);
	int *alle = INTEGER(alleles);
	for(i = 0; i < rows; i++)
	{
		// loop through all columns first and initialize
		for(j = 0; j < cols; j++) 
		{
			// maintain missing values
			if (ISNA(loc[i + j*rows]))
			{
				outmat[i + j*rows] = loc[i + j*rows]; 
				miss = 1;
			}
			else
			{
				outmat[i + j*rows] = 0.0; 
			}
		}
		if (miss == 1)
		{
			miss = 0;		 
		}
		else
		{
			// permute the alleles by adding the allele frequency
			for(p = 0; p < ploid; p++)
			{
				outmat[i + alle[count++]*rows] += *(afreq);
			}
		}
	}
	UNPROTECT(1);
	return Rout;
}

