#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP association_index_diploid(SEXP, SEXP, SEXP, SEXP);
extern SEXP association_index_haploid(SEXP, SEXP, SEXP);
extern SEXP bitwise_distance_diploid(SEXP, SEXP, SEXP, SEXP);
extern SEXP bitwise_distance_haploid(SEXP, SEXP, SEXP);
extern SEXP bruvo_distance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP expand_indices(SEXP, SEXP);
extern SEXP genotype_curve_internal(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_pgen_matrix_genind(SEXP, SEXP, SEXP, SEXP);
extern SEXP mlg_round_robin(SEXP);
extern SEXP msn_tied_edges(SEXP, SEXP, SEXP);
extern SEXP neighbor_clustering(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP omp_test();
extern SEXP pairdiffs(SEXP);
extern SEXP pairwise_covar(SEXP);
extern SEXP permute_shuff(SEXP, SEXP, SEXP);
extern SEXP permuto(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"association_index_diploid", (DL_FUNC) &association_index_diploid, 4},
    {"association_index_haploid", (DL_FUNC) &association_index_haploid, 3},
    {"bitwise_distance_diploid",  (DL_FUNC) &bitwise_distance_diploid,  4},
    {"bitwise_distance_haploid",  (DL_FUNC) &bitwise_distance_haploid,  3},
    {"bruvo_distance",            (DL_FUNC) &bruvo_distance,            6},
    {"expand_indices",            (DL_FUNC) &expand_indices,            2},
    {"genotype_curve_internal",   (DL_FUNC) &genotype_curve_internal,   4},
    {"get_pgen_matrix_genind",    (DL_FUNC) &get_pgen_matrix_genind,    4},
    {"mlg_round_robin",           (DL_FUNC) &mlg_round_robin,           1},
    {"msn_tied_edges",            (DL_FUNC) &msn_tied_edges,            3},
    {"neighbor_clustering",       (DL_FUNC) &neighbor_clustering,       5},
    {"omp_test",                  (DL_FUNC) &omp_test,                  0},
    {"pairdiffs",                 (DL_FUNC) &pairdiffs,                 1},
    {"pairwise_covar",            (DL_FUNC) &pairwise_covar,            1},
    {"permute_shuff",             (DL_FUNC) &permute_shuff,             3},
    {"permuto",                   (DL_FUNC) &permuto,                   1},
    {NULL, NULL, 0}
};

void R_init_poppr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
