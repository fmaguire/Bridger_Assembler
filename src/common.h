#ifndef COMMON_H
#define COMMON_H

/*
 *  common.h
 */

#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/typesizes.h>
#include <bits/types.h>


// general options
extern int g_kmer_length;
extern bool g_help;
extern bool g_debug;
extern bool g_double_stranded_mode;
extern bool g_is_paired_end;
extern std::string out_dir;
extern int g_fr_strand;

// assemble options
extern std::string reads_file;
extern std::string kmers_file;
extern int g_min_kmer_coverage;
extern float g_min_kmer_entropy;
extern int g_min_seed_coverage;
extern float g_min_seed_entropy;
extern int g_min_junction_coverage;
extern int g_min_average_coverage;
extern int g_min_anchor_length;
extern int g_min_reads_span_junction;
extern int g_min_reads_support_branch;
extern float g_min_ratio_non_error;
extern float g_min_ratio_branch;
extern float g_min_ratio_welds;
extern float g_min_ratio_in_out;
extern int g_min_exon_length;
extern int g_min_trunk_length;
extern int g_pair_gap_length;
extern int g_max_pair_gap_length;
extern int g_min_kmers_per_graph;
extern bool g_double_stranded_mode;
extern std::string rg_list;
extern int g_interval;

// path search options
extern int CPU;
extern std::string rg_file;
extern int g_min_transcript_length;
extern std::string output_filename;


#define OPT_MIN_KMER_COVERAGE		301
#define OPT_MIN_KMER_ENTROPY		302
#define OPT_MIN_SEED_COVERAGE		303
#define OPT_MIN_SEED_ENTROY		304
#define	OPT_MIN_JUNCTION_COVERAGE	305
#define OPT_MIN_RATIO_NON_ERROR		306
#define OPT_MIN_EXON_LENGTH		307
#define OPT_MIN_KMERS_PER_GRAPH		308
#define OPT_MIN_TRANSCRIPT_LENGTH	309
#define OPT_KMERS			311
#define OPT_DOUBLE_STRANDED_MODE	312
#define OPT_DEBUG			313
#define OPT_IS_PAIR_END			314
#define OPT_MAX_PAIR_GAP_LENGTH		315
#define OPT_PAIR_GAP_LENGTH		316
#define OPT_FR_STRAND			317

int mkpath(const char *s, mode_t mode);

#endif

