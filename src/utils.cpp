# include "utils.hpp"

namespace fs = std::filesystem;

int discrete_distn(const rundle& r) {
    std::vector<float> probs;
    probs.reserve(r.alns.size());
    for (const aln& a : r.alns) {
        assert (a.prob != -1);
        probs.push_back(a.prob);
    }

    static std::random_device seed;
    static std::mt19937 gen(seed());

    std::discrete_distribution<> dist(probs.begin(), probs.end());
    return dist(gen);
}

xt::xarray<float> assign_by_sampling(std::vector<rundle>& rundles,
                    std::vector<size_t>& multi_mappers, size_t tsize) {

    xt::xarray<float> new_rho = xt::zeros<float>({tsize});
    int num_threads = omp_get_max_threads();
    std::vector<xt::xarray<float>> thread_rhos(num_threads, xt::zeros<float>({tsize}));

    #pragma omp parallel
    {
        int pid = omp_get_thread_num();
        xt::xarray<float>& local_rho = thread_rhos[pid];
        # pragma omp for schedule(static)
        for (size_t mi = 0; mi < multi_mappers.size(); ++mi) {
            size_t i = multi_mappers[mi];
            rundle& r = rundles[i];
            int sel_ai = discrete_distn(r);
            for (size_t j = 0; j < r.alns.size(); ++j) {
                aln& a = r.alns[j];
                a.prob = (j == static_cast<size_t>(sel_ai)) ? 1.0f : 0.0f;
                local_rho[a.tid] += a.prob;
            }
        }

    }
    for (int t = 0; t < num_threads; ++t) {
        new_rho += thread_rhos[t];
    }
    return new_rho;
}

bool is_fastq(const std::string& fn) {
    size_t len = fn.size();
    return (len >= 6 && fn.compare(len - 6, 6, ".fastq") == 0) ||
           (len >= 3 && fn.compare(len - 3, 3, ".fq") == 0) ||
           (len >= 9 && fn.compare(len - 9, 9, ".fastq.gz") == 0) ||
           (len >= 6 && fn.compare(len - 6, 6, ".fq.gz") == 0);
};

bool is_bam(const std::string& fn) {
    size_t len = fn.size();
    return (len >= 6 && fn.compare(len - 4, 4, ".bam") == 0) ||
           (len >= 3 && fn.compare(len - 7, 7, ".bam.gz") == 0);
};

std::string join_path(const std::string& out_dir, const std::string& fn) {
    fs::path dir(out_dir);
    fs::path file(fn);
    return (dir / file).string();
}

void write_abundances(const std::string& fn, \
                    const xt::xarray<float>& final_rho, \
                    char** tnames, size_t tsize) {
    
    std::vector<std::string> lines(tsize);

    # pragma omp parallel for
    for (size_t i = 0; i < tsize; ++i) {
        std::ostringstream ss;
        ss << tnames[i] << "," << final_rho[i];
        lines[i] = ss.str();
    }

    std::ofstream outfile(fn);
    for (const auto& line : lines) {
        outfile << line << '\n';
    }

    outfile.close();
}

void write_assignments(const std::string& fn, std::vector<rundle>& rundles, \
                    char** tnames, size_t qsize) {

    std::vector<std::string> lines(qsize);
    std::string one_str = "1.0";

    # pragma omp parallel for
    for (size_t i = 0; i < qsize; ++i) {
        const rundle& r = rundles[i];
        std::ostringstream ss;
        ss << r.qname;

        for (const aln& a : r.alns) {
            if (a.prob == 0) {
                continue;
            } else {
                ss << '\t' << tnames[a.tid] << ',';
                if (a.prob == -1) {
                    ss << one_str; // uniquely aligned
                } else {
                    ss << a.prob;
                }
            }
        }

        lines[i] = ss.str();
    }

    std::ofstream outfile(fn);
    for (const auto& line : lines) {
        outfile << line << '\n';
    }

    outfile.close();
}

void drop_low_prob(std::vector<rundle>& rundles, \
                std::vector<size_t>& multi_mappers, \
                const float min_pk, const float base_w) {
    # pragma omp for schedule(static)
    for (size_t mi = 0; mi < multi_mappers.size(); ++mi) {
        size_t i = multi_mappers[mi];
        rundle& r = rundles[i];
        float min_prob = (1.0f + min_pk) / r.alns.size();
        for (aln& a : r.alns) {
            if (a.score == r.max_score) {
                continue;
            } else {
                if (a.prob < min_prob) {
                    a.weight = base_w;
                }
            }
        }
    }
}

void zero_out_rho(xt::xarray<float>& rho, float min_rho) {
    rho = xt::where(rho < min_rho, 0.0f, rho);
}