/*
 *  common.h
 */
#include "common.h"
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <bits/types.h>
//#include <bits/typesizes.h>
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <libgen.h>

// general
int g_kmer_length = 25;
bool g_help = false;
bool g_debug = false;
bool g_is_paired_end = false;
int g_fr_strand = 1;

// assemble
int g_min_kmer_coverage = 1;
float g_min_kmer_entropy = 0.0f;
int g_min_seed_coverage = 2;
float g_min_seed_entropy = 1.5f;
int g_min_junction_coverage = 2;
int g_min_average_coverage = 2;
int g_min_reads_span_junction = 2;
int g_min_reads_support_branch = 2;
int g_min_anchor_length = 21;
float g_min_ratio_non_error = 0.05f;
float g_min_ratio_branch = 0.05f;
float g_min_ratio_welds = 0.04f;
float g_min_ratio_in_out = 0.02f;
int g_min_exon_length = 80;
int g_min_trunk_length = 200;
int g_min_kmers_per_graph = 276;
int g_pair_gap_length = 200;
int g_max_pair_gap_length = 500;
std::string reads_file = "";
std::string kmers_file = "";
std::string out_dir = "./RawGraphs/";
bool g_double_stranded_mode = false;
std::string rg_list = "raw_graph.list";
int g_interval = 20;

// path search
int CPU = 6;
std::string rg_file = "";
int g_min_transcript_length = 200;
std::string output_filename = "transcripts.fasta";


/* Function with behaviour like `mkdir -p'  */
/* found at: http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/ */
/*
int mkpath(const char *s, mode_t mode) {
    char *q, *r = NULL, *path = NULL, *up = NULL;
    int rv;
    rv = -1;
    if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
        return (0);
    if ((path = strdup(s)) == NULL)
        exit(1);
    if ((q = strdup(s)) == NULL)
        exit(1);
    if ((r = dirname(q)) == NULL)
        goto out;
    if ((up = strdup(r)) == NULL)
        exit(1);
    if ((mkpath(up, mode) == -1) && (errno != EEXIST))
        goto out;
    if ((mkdir(path, mode) == -1) && (errno != EEXIST))
        rv = -1;
    else
        rv = 0;
out:
    if (up != NULL)
        free(up);
    free(q);
    free(path);
    return (rv);
}
*/
