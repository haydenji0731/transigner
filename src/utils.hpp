# ifndef UTILS_H
# define UTILS_H

# include <xtensor/containers/xarray.hpp>
# include <xtensor/io/xio.hpp>
# include <xtensor/io/xio.hpp>
# include <xtensor/views/xview.hpp>
# include <xtensor/core/xmath.hpp>
# include <fstream>
# include <algorithm>
# include <omp.h>
# include <stdlib.h>
# include <iostream>
# include <string>
# include <unordered_set>
# include <unordered_map>
# include <vector>
# include <string>
# include "htslib/sam.h"
# include <chrono>
# include <cstdio>
# include <cstdlib>
# include <cmath>
# include <filesystem>
# include <fstream>
# include <random>

# define RED_TEXT(msg) "\033[1;31m" msg "\033[0m"
# define GREEN_TEXT(msg) "\033[1;32m" msg "\033[0m"

const int MAX_N_READS = 10000000; // ~10 mio reads
const int MAX_N_ALNS = 100;

struct aln {
    int32_t tid;
    hts_pos_t start, end;
    mutable float score;
    mutable float weight = -1;
    mutable float prob;
};

struct rundle {
    char* qname;
    std::vector<aln> alns;
    float max_score;
};

void write_assignments(const std::string& fn, std::vector<rundle>& rundles, \
                    char** tnames, size_t qsize);
void write_abundances(const std::string& fn, const xt::xarray<float>& final_rho, \
                    char** tnames, size_t tsize);

// TODO
int discrete_distn(const rundle& r);
xt::xarray<float> assign_by_sampling(std::vector<rundle>& rundles,
                    std::vector<size_t>& multi_mappers, size_t tsize);

void drop_low_prob(std::vector<rundle>& rundles,
                std::vector<size_t>& multi_mappers,
                const float min_pk, const float base_w);
void zero_out_rho(xt::xarray<float>& rho, float thres);

bool is_fastq(const std::string& fn);
bool is_bam(const std::string& fn);
std::string join_path(const std::string& out_dir, const std::string& fn);

# endif