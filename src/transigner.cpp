# include "utils.hpp"
# include <cxxopts.hpp>

# define VERSION "1.2.0"

cxxopts::ParseResult parse_args(int argc, char* argv[]) {
    cxxopts::Options options("transigner", "TranSigner version " VERSION);

    // bool values don't display the default_value automatically
    options.add_options()
        ("o,out-dir", "output directory", cxxopts::value<std::string>()->default_value("."))
        ("t,tx", "fasta file containing transcript sequences", cxxopts::value<std::string>())
        ("a,max-aln", "maximum number of alignments per read", cxxopts::value<int>()->default_value("181"))
        ("n,max-iter", "number of max EM iterations", cxxopts::value<int>()->default_value("1000"))
        ("e,epsilon", "EM convergence threshold ε", cxxopts::value<float>()->default_value("10.0"))
        ("c,constant", "τ constant in exponential decay fn", cxxopts::value<float>()->default_value("5.0"))
        ("sample", "assign by sampling (default: false)", cxxopts::value<bool>()->default_value("false"))
        ("use-psw", "use positional weights (default: false)", cxxopts::value<bool>()->default_value("false"))
        ("r,min-rho", "zero out txes with read counts below this threshold", cxxopts::value<float>()->default_value("0.1"))
        ("keep-low", "keep low-prob alignments (default: false)", cxxopts::value<bool>()->default_value("false"))
        ("k,min-pk", "constant to computed min alignment prob", cxxopts::value<float>()->default_value("0.1"))
        ("w,base-w", "base weight for an alignment", cxxopts::value<float>()->default_value("0.0"))
        ("d,data-type", "long read RNA-seq data type", cxxopts::value<std::string>()->default_value("ont"))
        ("p,threads", "number of threads", cxxopts::value<int>()->default_value("1"))
        ("v,verbose", "enable verbose output for more detailed logging (default: false)", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "usage")
        ("V,version", "Print version information and exit")
        ("input", "input filename", cxxopts::value<std::string>())
    ;
    options.parse_positional({"input"});

    auto res = options.parse(argc, argv);

    if (res.count("version")) {
        std::cout << "transigner " << GREEN_TEXT("v") << GREEN_TEXT(VERSION) << std::endl;
        std::exit(0);
    }

    if (res.count("help")) {
        std::cout << options.help() << std::endl;
        std::exit(EXIT_SUCCESS);
    }
    return res;
}

void align(const std::string reads_file, const std::string txome_file, \
    const std::string out_file, const std::string dtype, size_t threads, \
    size_t num_alignments) {
    
    std::string preset;

    if (dtype == "ont") {
        preset = "map-ont";
    } else if (dtype == "pacbio") {
        preset = "map-pb";
    } else {
        std::cerr << RED_TEXT("error: ") << "unknown data type '" << dtype << "'\n";
    } // TODO: consider adding isoseq specific preset

    std::string cmd = "minimap2 -ax " + preset + " -N " + std::to_string(num_alignments) + 
                    " -t " + std::to_string(threads) + " " + txome_file + " " +
                    reads_file + " | samtools view -b -o " + out_file +
                    " -@ " + std::to_string(threads);
    std::cout << cmd << std::endl;
    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << RED_TEXT("error: ") << "alignment command failed with exit code " << ret << "\n";
        std::exit(EXIT_FAILURE);
    }
}

std::vector<rundle> load_alignments(htsFile* hts_file, sam_hdr_t* hdr, bam1_t* b_next, \
    std::unordered_set<int32_t>& aligned_tids){
    std::vector<rundle> rundles;
    rundles.reserve(MAX_N_READS);

    std::vector<aln> alns;
    alns.reserve(MAX_N_ALNS);
    int max_score = -1;
    char* prev_name = nullptr;

    while (sam_read1(hts_file, hdr, b_next) >= 0) {
        if (((b_next->core.flag & BAM_FUNMAP) != 0) || \
            ((b_next->core.flag & BAM_FSUPPLEMENTARY) != 0)) {
            continue;
        }
        const char* curr_name = bam_get_qname(b_next);
        if (prev_name == nullptr) {
            prev_name = strdup(curr_name);
        } else if (strcmp(prev_name, curr_name) != 0) {
            rundle r;
            r.qname = strdup(prev_name);
            r.alns = alns;
            r.max_score = max_score;
            rundles.push_back(std::move(r));

            free(prev_name);
            prev_name = strdup(curr_name);
            alns.clear();
            max_score = -1;
        }
        aln a;
        a.tid = b_next->core.tid;
        aligned_tids.insert(a.tid);
        a.start = b_next->core.pos;
        a.end = bam_endpos(b_next);
        a.score = bam_aux2i(bam_aux_get(b_next, "AS"));
        if (a.score > max_score) {
            max_score = a.score;
        }
        alns.push_back(std::move(a));
    }

    if (alns.size() >= 1) {
        rundle r;
        r.qname = strdup(prev_name);
        r.alns = alns;
        r.max_score = max_score;
        rundles.push_back(std::move(r));
    }
    return rundles;
}

std::vector<size_t> compute_weights_psw(std::vector<rundle>& rundles, size_t tsize, \
    size_t qsize, uint32_t* tlens, std::vector<xt::xarray<float>>& pos_ws, \
    std::vector<xt::xarray<float>>& score_ws, xt::xarray<float>& uniq_rho) {
    
    xt::xarray<float> tcov = xt::zeros<float>({tsize});
    std::vector<size_t> multi_mappers;
    multi_mappers.reserve(qsize);

    for (size_t i = 0; i < qsize; ++i) {
        auto& rundle = rundles[i];
        if (rundle.alns.size() == 1) { // uniquely aligned
            uniq_rho[rundle.alns[0].tid]++;
            rundle.alns[0].prob = -1;
            continue;
        }
        multi_mappers.push_back(i);

        xt::xarray<float> score_data = xt::empty<float>({rundle.alns.size()});
        for (size_t j = 0; j < rundle.alns.size(); ++j) {
            auto& a = rundle.alns[j];
            int32_t tid = a.tid;
            hts_pos_t start = a.start;
            hts_pos_t end = a.end;
            tcov[tid] += end - start;
            if (pos_ws[tid].size() == 1) {
                pos_ws[tid] = xt::zeros<float>({tlens[tid]});
            }
            xt::view(pos_ws[tid], xt::range(start, end)) += 1;
            score_data[j] = a.score; 
        }
        score_data = (score_data - rundle.max_score) / 5.0f; // TODO: expose this param
        score_ws[i] = xt::exp(score_data);
    }
    for (size_t k = 0; k < tsize; ++k) {
        if (tcov[k] == 0) {
            continue; // no alignment
        }
        auto& pos_ws_k = pos_ws[k];
        float norm_tcov = tcov[k] / static_cast<float>(tlens[k]);
        
        pos_ws_k = norm_tcov - pos_ws_k;
        pos_ws_k = xt::eval(pos_ws_k);
    
        float min_val = std::abs(xt::amin(pos_ws_k)());
        pos_ws_k += min_val;
    
        float z = xt::sum(pos_ws_k)();
        if (z > 0) {
            pos_ws_k = xt::cumsum(pos_ws_k / z, 0);
        } else {
            float inv_len = 1.0f / tlens[k];
            pos_ws_k = xt::cumsum(xt::ones<float>({tlens[k]}) * inv_len, 0);
        }
    }
    return multi_mappers;
}

std::vector<size_t> compute_weights(std::vector<rundle>& rundles, size_t qsize, \
    std::vector<xt::xarray<float>>& score_ws, xt::xarray<float>& uniq_rho) {
    
    std::vector<size_t> multi_mappers;
    multi_mappers.reserve(qsize);

    for (size_t i = 0; i < qsize; ++i) {
        auto& rundle = rundles[i];
        if (rundle.alns.size() == 1) { // uniquely aligned
            uniq_rho[rundle.alns[0].tid]++;
            rundle.alns[0].prob = -1;
            continue;
        }
        multi_mappers.push_back(i);

        xt::xarray<float> score_data = xt::empty<float>({rundle.alns.size()});
        for (size_t j = 0; j < rundle.alns.size(); ++j) {
            auto& a = rundle.alns[j];
            score_data[j] = a.score; 
        }
        score_data = (score_data - rundle.max_score) / 5.0f; // TODO: expose this
        score_ws[i] = xt::exp(score_data);
    }
    return multi_mappers;
}

xt::xarray<float> do_em_init_psw(std::vector<size_t>& multi_mappers, std::vector<rundle>& rundles, \
    std::vector<xt::xarray<float>>& pos_ws, std::vector<xt::xarray<float>>& score_ws, \
    xt::xarray<float>& rho, size_t tsize) {

    xt::xarray<float> new_rho = xt::zeros<float>({tsize});
    int num_threads = omp_get_max_threads();
    std::vector<xt::xarray<float>> thread_rhos(num_threads, xt::zeros<float>({tsize}));

    # pragma omp parallel
    {
        int pid = omp_get_thread_num();
        xt::xarray<float>& local_rho = thread_rhos[pid];

        # pragma omp for schedule(static)
        for (size_t mi = 0; mi < multi_mappers.size(); ++mi) {
            size_t i = multi_mappers[mi];
            rundle& r = rundles[i];

            float z = 0.0f;
            for (const aln& a : r.alns) {
                float weight;
                if (a.weight == -1) {
                    weight = score_ws[i][&a - &r.alns[0]] * \
                        (pos_ws[a.tid][a.end - 1] - pos_ws[a.tid][a.start]);
                } else {
                    weight = a.weight;
                }
                z += rho[a.tid] * weight;
            }

            if (z > 0.0f) {
                for (aln& a : r.alns) {
                    float weight;
                    if (a.weight == -1) {
                        weight = score_ws[i][&a - &r.alns[0]] * \
                            (pos_ws[a.tid][a.end - 1] - pos_ws[a.tid][a.start]);
                        a.weight = weight;
                    } else {
                        weight = a.weight;
                    }
                    a.prob = (rho[a.tid] * weight) / z;
                    local_rho[a.tid] += a.prob;
                }
            }
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        new_rho += thread_rhos[t];
    }
    return new_rho;
}

xt::xarray<float> do_em_init(std::vector<size_t>& multi_mappers, std::vector<rundle>& rundles, \
    std::vector<xt::xarray<float>>& score_ws, xt::xarray<float>& rho, size_t tsize) {

    xt::xarray<float> new_rho = xt::zeros<float>({tsize});
    int num_threads = omp_get_max_threads();
    std::vector<xt::xarray<float>> thread_rhos(num_threads, xt::zeros<float>({tsize}));

    # pragma omp parallel
    {
        int pid = omp_get_thread_num();
        xt::xarray<float>& local_rho = thread_rhos[pid];

        # pragma omp for schedule(static)
        for (size_t mi = 0; mi < multi_mappers.size(); ++mi) {
            size_t i = multi_mappers[mi];
            rundle& r = rundles[i];

            float z = 0.0f;
            for (const aln& a : r.alns) {
                float weight;
                if (a.weight == -1) {
                    weight = score_ws[i][&a - &r.alns[0]];
                } else {
                    weight = a.weight;
                }
                z += rho[a.tid] * weight;
            }

            if (z > 0.0f) {
                for (aln& a : r.alns) {
                    float weight;
                    if (a.weight == -1) {
                        weight = score_ws[i][&a - &r.alns[0]];
                        a.weight = weight;
                    } else {
                        weight = a.weight;
                    }
                    a.prob = (rho[a.tid] * weight) / z;
                    local_rho[a.tid] += a.prob;
                }
            }
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        new_rho += thread_rhos[t];
    }
    return new_rho;
}

xt::xarray<float> do_em(std::vector<size_t>& multi_mappers, std::vector<rundle>& rundles, \
    xt::xarray<float>& rho, size_t tsize) {

    xt::xarray<float> new_rho = xt::zeros<float>({tsize});
    int num_threads = omp_get_max_threads();
    std::vector<xt::xarray<float>> thread_rhos(num_threads, xt::zeros<float>({tsize}));
    
    #pragma omp parallel
    {
        int pid = omp_get_thread_num();
        xt::xarray<float>& local_rho = thread_rhos[pid];

        #pragma omp for schedule(static)
        for (size_t mi = 0; mi < multi_mappers.size(); ++mi) {
            size_t i = multi_mappers[mi];
            rundle& r = rundles[i];

            float z = 0.0f;
            for (const aln& a : r.alns) {
                z += rho[a.tid] * a.weight;
            }

            if (z > 0.0f) {
                for (aln& a : r.alns) {
                    a.prob = (rho[a.tid] * a.weight) / z;
                    local_rho[a.tid] += a.prob;
                }
            }
        }
    }
    for (int t = 0; t < num_threads; ++t) {
        new_rho += thread_rhos[t];
    }
    return new_rho;
}

int main(int argc, char* argv[]) {

    auto res = parse_args(argc, argv);

    std::string in_fn = res["input"].as<std::string>();
    int num_threads = res["threads"].as<int>();
    int num_alignments = res["max-aln"].as<int>();
    std::string out_dir = res["out-dir"].as<std::string>();

    if (is_fastq(in_fn)) {
        if (!res.count("tx")) {
            std::cerr << RED_TEXT("error: ") << "missing transcriptome file" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const std::string& txome_fn = res["tx"].as<std::string>();
        const std::string& dtype = res["data-type"].as<std::string>();
        std::string out_fn = join_path(out_dir, "alignments.bam");

        align(in_fn, txome_fn, out_fn, dtype, num_threads, num_alignments);
        in_fn = out_fn;
        
    } else if(is_bam(in_fn)) {
        std::cout << GREEN_TEXT("status: ") << "skipping alignment" << std::endl;
    } else {
        std::cerr << RED_TEXT("error: ") << "unknown input file type" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    omp_set_num_threads(num_threads);

    htsFile* hts_file = hts_open(in_fn.c_str(), "r");
    sam_hdr_t* hdr = sam_hdr_read(hts_file);
    bam1_t* b_next = bam_init1();

    size_t tsize = hdr->n_targets;
    uint32_t* tlens = hdr->target_len;
    char** tnames = hdr->target_name;

    std::unordered_set<int32_t> aligned_tids;
    aligned_tids.reserve(tsize);

    auto rundles = load_alignments(hts_file, hdr, b_next, aligned_tids);
    
    size_t qsize = rundles.size();

    const bool use_psw = res["use-psw"].as<bool>();
    if (use_psw) {
        std::cout << GREEN_TEXT("status: ") << "using position-specific weights" << std::endl;
    }

    std::vector<xt::xarray<float>> score_ws(qsize);
    std::vector<xt::xarray<float>> pos_ws;
    xt::xarray<float> uniq_rho = xt::zeros<float>({tsize});
    std::vector<size_t> multi_mappers;

    if (use_psw) {
        pos_ws.resize(tsize);
        multi_mappers = compute_weights_psw(rundles, tsize, qsize, \
            tlens, pos_ws, score_ws, uniq_rho);
    } else {
        multi_mappers = compute_weights(rundles, qsize, score_ws, uniq_rho);
    }

    xt::xarray<float> rho = xt::zeros<float>({tsize});
    // TODO: consider initializing just the aligned tids; shouldn't have much impact on EM results
    rho.fill(static_cast<float>(qsize) / static_cast<float>(aligned_tids.size()));

    int ctr = 0;
    xt::xarray<float> new_rho;
    float delta = 0.0f;
    const float epsilon = res["epsilon"].as<float>();
    const int max_iter = res["max-iter"].as<int>();
    const float min_rho = res["min-rho"].as<float>();
    bool keep_low_prob = res["keep-low"].as<bool>();
    if (use_psw && !keep_low_prob) {
        keep_low_prob = true;
    }
    const float min_pk = res["min-pk"].as<float>();
    const float base_w = res["base-w"].as<float>();
    
    std::cout << GREEN_TEXT("status: ") << "em started" << std::endl;

    while (ctr < max_iter) {

        # ifdef DEBUG
        auto start = std::chrono::high_resolution_clock::now();
        # endif

        if (ctr == 0) {
            if (use_psw) {
                new_rho = do_em_init_psw(multi_mappers, rundles, pos_ws, score_ws, rho, tsize);
            } else {
                new_rho = do_em_init(multi_mappers, rundles, score_ws, rho, tsize);
            }
        } else {
            new_rho = do_em(multi_mappers, rundles, rho, tsize);
        }

        if (ctr == 0 && !keep_low_prob) {
            drop_low_prob(rundles, multi_mappers, min_pk, base_w);
            new_rho = do_em(multi_mappers, rundles, rho, tsize);
        }

        new_rho += uniq_rho;
        
        if (min_rho > 0.0) {
            zero_out_rho(new_rho, min_rho);
        }

        delta = xt::sum(xt::abs(new_rho - rho))();

        if (delta < epsilon) {
            std::cout << GREEN_TEXT("status: ") << "converged @ it " << ctr << std::endl;
            break;
        }

        rho = std::move(new_rho);

        # ifdef DEBUG
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        printf("\tit %zu: %.6f seconds\n", ctr, elapsed.count());
        # endif

        ctr++;
    }

    if (ctr == max_iter) {
        std::cout << GREEN_TEXT("status: ") << "max iteration (" << max_iter << ") reached" << std::endl;
    }

    const bool do_sample = res["sample"].as<bool>();

    if (do_sample) {
        new_rho = assign_by_sampling(rundles, multi_mappers, tsize);
        new_rho += uniq_rho;

        # ifdef DEBUG
        delta = xt::sum(xt::abs(new_rho - rho))();
        printf("\tdelta after sampling: %f\n", delta);
        # endif

        rho = std::move(new_rho);
    }

    write_abundances(join_path(out_dir, "abundances.csv"), rho, tnames, tsize);
    write_assignments(join_path(out_dir, "assignments.tsv"), rundles, tnames, qsize);

    bam_destroy1(b_next);
    hts_close(hts_file);
    sam_hdr_destroy(hdr);

    return 0;
}