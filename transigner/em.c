#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "xxhash.h"
#ifdef _OPENMP
    #include <omp.h>
#endif

#include <time.h>

#define HMAP_SIZE 150000
#define MAX_LINE_LENGTH 1024
#define MAX_FIELDS 100
#define PATH_SEPARATOR '/'

typedef const unsigned int cuint;
typedef const unsigned long culn;
typedef const unsigned long long culnln;
typedef unsigned int uint;

typedef int bool;
#define true 1
#define false 0
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct opt {
    char *flag;
    char *desc;
    char *value;
    bool requires_val;
} opt_t;

typedef struct kv_pair {
    char *key;
    void *value;
    struct kv_pair *next;
} kv_pair_t;

typedef struct {
    kv_pair_t **buckets;
} hmap_t;

unsigned long long hash_fn(const char *key) {
    return XXH64(key, strlen(key), 0); // use 0 as the seed
}

hmap_t *hmap_init() {
    hmap_t *hmap = malloc(sizeof(hmap_t));
    hmap->buckets = malloc(sizeof(kv_pair_t *) * HMAP_SIZE);
    for (uint bi = 0; bi < HMAP_SIZE; bi++) {
        hmap->buckets[bi] = NULL;
    }
    return hmap;
}

void hmap_insert(hmap_t *hmap, const char *key, void *value) {
    unsigned long long hash = hash_fn(key);
    uint index = hash % HMAP_SIZE;
    kv_pair_t *new_pair = malloc(sizeof(kv_pair_t));
    new_pair->key = strdup(key);
    new_pair->value = value;
    new_pair->next = hmap->buckets[index]; // resolve hashing collision by chaining
    hmap->buckets[index] = new_pair;
}

void *hmap_lookup(hmap_t *hmap, const char *key) {
    unsigned long long hash = hash_fn(key);
    int index = hash % HMAP_SIZE;
    kv_pair_t *curr = hmap->buckets[index];
    while (curr != NULL) {
        if (strcmp(curr->key, key) == 0) {
            return curr->value;
        }
        curr = curr->next;
    }
    return NULL;
}

void hmap_free(hmap_t *hmap) {
    for (int bi = 0; bi < HMAP_SIZE; bi++) {
        kv_pair_t *curr = hmap->buckets[bi];
        while (curr != NULL) {
            kv_pair_t *temp = curr;
            curr = curr->next;
            free(temp->key);
            free(temp);
        }
    }
    free(hmap->buckets);
    free(hmap);
}

void hmap_print(const hmap_t *hmap) {
    for (int bi = 0; bi < HMAP_SIZE; bi++) {
        kv_pair_t *curr = hmap->buckets[bi];
        if (curr != NULL) {
            printf("Bucket %d:\n", bi);
            while (curr != NULL) {
                printf("\tKey: %s, Value: %d\n", curr->key, *(int *)curr->value);
                curr = curr->next;
            }
        }
    }
}

typedef struct tinfo_t {
    char * tname;
    double score;
    double rfrac;
} tinfo_t;

typedef struct qinfo_t {
    char * qname;
    tinfo_t ** tinfo_vec;
    uint n;
} qinfo_t;

void parse_csv_line(char *line, char fields[MAX_FIELDS][MAX_LINE_LENGTH], int *field_count) {
    char *token = strtok(line, ",");
    *field_count = 0;
    while (token != NULL) {
        strncpy(fields[*field_count], token, MAX_LINE_LENGTH);
        (*field_count)++;
        token = strtok(NULL, ",");
    }
}

qinfo_t** load_qt_info(const char* ifn, cuint qsize, cuint tsize) {
    FILE *file = fopen(ifn, "r");
    char line[MAX_LINE_LENGTH];
    char fields[MAX_FIELDS][MAX_LINE_LENGTH];
    int field_count;
    bool is_first = true;
    uint i = 0;
    uint qi = 0;
    uint ti = 0;
    qinfo_t *qinfo = (qinfo_t*) malloc(sizeof(qinfo_t));
    tinfo_t *tinfo = (tinfo_t*) malloc(sizeof(tinfo_t));
    qinfo_t **qinfo_vec = (qinfo_t**) malloc(qsize * sizeof(qinfo_t*));
    fgets(line, sizeof(line), file); // skip the first line
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = '\0';
        parse_csv_line(line, fields, &field_count);
        if (is_first) {
            is_first = false;
            qinfo->n = 1;
            qinfo->qname = (char*) malloc(strlen(fields[0]) + 1);
            strcpy(qinfo->qname, fields[0]);
            qinfo->tinfo_vec = (tinfo_t**) malloc(tsize * sizeof(tinfo_t*));
            tinfo->tname = (char*) malloc(strlen(fields[1]) + 1);
            strcpy(tinfo->tname, fields[1]);
            tinfo->score = strtod(fields[2], NULL);
            tinfo->rfrac = 0.0;
            qinfo->tinfo_vec[ti] = tinfo;
            ti++;
        } else if (strcmp(qinfo->qname, fields[0]) != 0) {
            qinfo_vec[qi] = qinfo;
            qi++;
            qinfo = (qinfo_t*) malloc(sizeof(qinfo_t));
            ti = 0;
            qinfo->n = 1;
            qinfo->qname = (char*) malloc(strlen(fields[0]) + 1);
            strcpy(qinfo->qname, fields[0]);
            qinfo->tinfo_vec = (tinfo_t**) malloc(tsize * sizeof(tinfo_t*));
            tinfo = (tinfo_t*) malloc(sizeof(tinfo_t));
            tinfo->tname = (char*) malloc(strlen(fields[1]) + 1);
            strcpy(tinfo->tname, fields[1]);
            tinfo->score = strtod(fields[2], NULL);
            tinfo->rfrac = 0.0;
            qinfo->tinfo_vec[ti] = tinfo;
            ti++;
        } else {
            tinfo = (tinfo_t*) malloc(sizeof(tinfo_t));
            tinfo->tname = (char*) malloc(strlen(fields[1]) + 1);
            strcpy(tinfo->tname, fields[1]);
            tinfo->score = strtod(fields[2], NULL);
            tinfo->rfrac = 0.0;
            qinfo->tinfo_vec[ti] = tinfo;
            ti++;
            qinfo->n++;
        }
    }
    qinfo_vec[qi] = qinfo;
    qi++;
    // sanity check
    assert(qi == qsize);
    return qinfo_vec;
}

void free_qt_info(qinfo_t** qinfo_vec, cuint qsize, cuint tsize) {
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        if (qinfo != NULL) {
            if (qinfo->qname != NULL) {
                free(qinfo->qname);
            }
            if (qinfo->tinfo_vec != NULL) {
                for (uint ti = 0; ti < qinfo->n; ti++) {
                    tinfo_t *tinfo = qinfo->tinfo_vec[ti];
                    if (tinfo != NULL) {
                        if (tinfo->tname != NULL) {
                            free(tinfo->tname);
                        }
                    }
                    free(tinfo);
                }
            }
            free(qinfo->tinfo_vec);
        }
        free(qinfo);
    }
    free(qinfo_vec);
}

hmap_t * load_tmap(const char* ifn, char ** tnames) {
    FILE *file = fopen(ifn, "r");
    char line[MAX_LINE_LENGTH];
    char fields[MAX_FIELDS][MAX_LINE_LENGTH];
    int field_count;
    hmap_t *tmap = hmap_init();
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = '\0';
        parse_csv_line(line, fields, &field_count);
        int *value = malloc(sizeof(int));
        *value = (int) strtol(fields[1], NULL, 10);
        tnames[*value] = strdup(fields[0]);
        hmap_insert(tmap, fields[0], value);
    }
    fclose(file);
    return tmap;
}

void print_tnames(const char ** tnames, cuint tsize) {
    for (uint ti = 0; ti < tsize; ti++) {
        printf("%s\n", tnames[ti]);
    }
}

void print_qinfo_basic(const qinfo_t ** qinfo_vec, cuint qsize) {
    for (uint qi = 0; qi < qsize; qi++) {
        const qinfo_t *qinfo = qinfo_vec[qi];
        if (qinfo != NULL) {
            printf("qinfo[%d].n = %d\t", qi, qinfo->n);
            printf("qinfo[%d].qname = %s\n", qi, qinfo->qname);
            printf("qinfo[%d].score = %f\n", qi, qinfo->tinfo_vec[0]->score);
        } else {
            printf("qinfo[%d] is NULL\n", qi);
        }
    }
}

// TODO: observe if tsize has an impact?
double * init(cuint qsize, cuint tsize) {
    double * rho = (double*) malloc(tsize * sizeof(double));
    for (uint ti = 0; ti < tsize; ti++) {
        rho[ti] = qsize / tsize;
    }
    return rho;
}

void step_e(qinfo_t **qinfo_vec, double *rho, hmap_t *tmap, cuint qsize, bool is_naive) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        // printf("%s\n", qinfo->qname);
        if (qinfo == NULL) {
            printf("qinfo[%d] is NULL\n", qi);
            exit(-1); // an error
        }
        double z = 0.0;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            int *rho_i = (int *) hmap_lookup(tmap, tinfo->tname);
            if (rho_i == NULL) {
                printf("rho index is NULL\n");
                exit(-1); // an error
            }
            double rfrac;
            if (is_naive) {
                rfrac = rho[*rho_i];
            } else {
                rfrac = rho[*rho_i] * tinfo->score;
            }
            // assert(rho[*rho_i] > 0);
            tinfo->rfrac = rfrac;
            // printf("%s\t%.20lf\n", tinfo->tname, rfrac);
            z += rfrac;
        } // for now, allow no assignment
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            if (z == 0) {
                tinfo->rfrac = 0; // TODO: investigate them further
            } else {
                tinfo->rfrac /= z;
            }
            // if (isnan(tinfo->rfrac)) {
            //     printf("%s\t%.20lf\t%.20lf\n", tinfo->tname, tinfo->rfrac, z);
            // }
            // printf("%s\t%.20lf\n", tinfo->tname, tinfo->rfrac);
        }

    }
}

void step_m(qinfo_t **qinfo_vec, double *rho, hmap_t *tmap, cuint qsize) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            int *rho_i = (int *) hmap_lookup(tmap, tinfo->tname);
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            rho[*rho_i] += tinfo->rfrac;
        }
    }
}

void write_results(const qinfo_t **qinfo_vec, cuint qsize, cuint tsize, \
                double *rho, char * asgn_fn, char * unasgn_fn, char * abnd_fn, \
                char ** tnames) {
    FILE *asgn_file = fopen(asgn_fn, "w");
    FILE *unasgn_file = fopen(unasgn_fn, "w");
    for (uint qi = 0; qi < qsize; qi++) {
        const qinfo_t *qinfo = qinfo_vec[qi];
        double z = 0.0;
        bool first = true;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            z += tinfo->rfrac;
            if (tinfo->rfrac > 0.0) {
                if (first) {
                    fprintf(asgn_file, "%s,%s:%f", qinfo->qname, tinfo->tname, tinfo->rfrac);
                    first = false;
                } else {
                    fprintf(asgn_file, ",%s:%f", tinfo->tname, tinfo->rfrac);
                }
            }
        }
        if (z == 0.0) {
            fprintf(unasgn_file, "%s\n", qinfo->qname);
        } else {
            fprintf(asgn_file, "\n");
        }
    }
    fclose(asgn_file);
    fclose(unasgn_file);

    FILE *abnd_file = fopen(abnd_fn, "w");
    for (uint ti = 0; ti < tsize; ti++) {
        char * tname = tnames[ti];
        fprintf(abnd_file, "%s,%f\n", tname, rho[ti]);
    }
    fclose(abnd_file);
}

void print_help_msg(opt_t *cmd_opts, int n_opts) {
    printf("Usage: ./em [options]\n\n");
    printf("Options:\n");
    for (int i = 0; i < n_opts; i++) {
        printf("  %-15s %s\n", cmd_opts[i].flag, cmd_opts[i].desc);
    }
}

int is_int(const char *s) {
    char *endptr;
    strtol(s, &endptr, 10);
    return (*endptr == '\0');
}

int is_float(const char *s) {
    if (s == NULL || *s == '\0') {
        return 0;
    }
    char *endptr;
    strtod(s, &endptr);
    return (*endptr == '\0');
}

int discrete_distn(double * ws, int n, double z) {
    double rv = (double) rand() / (double) RAND_MAX;
    double running;
    for (int i = 0; i < n; i++) {
        running += ws[i];
        if (rv < running) {
            return i;
        }
    }
    return -1;
}

int parse_args(int argc, char *argv[], opt_t *cmd_opts, int n_opts) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help_msg(cmd_opts, n_opts);
            return 1;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--scores") == 0) {
            if (i + 1 < argc) {
                cmd_opts[1].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --scores option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--tmap") == 0) {
            if (i + 1 < argc) {
                cmd_opts[2].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --tmap option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-qs") == 0 || strcmp(argv[i], "--qsize") == 0) {
            if (i + 1 < argc) {
                if (is_int(argv[i + 1])) {
                    cmd_opts[3].value = argv[i + 1];
                    i++;
                } else {
                    fprintf(stderr, "Error: --qsize requires an integer value, but got '%s'.\n", argv[i + 1]);
                    return -1;
                }
            } else {
                fprintf(stderr, "Error: --qsize option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-ts") == 0 || strcmp(argv[i], "--tsize") == 0) {
            if (i + 1 < argc) {
                if (is_int(argv[i + 1])) {
                    cmd_opts[4].value = argv[i + 1];
                    i++;
                } else {
                    fprintf(stderr, "Error: --tsize requires an integer value, but got '%s'.\n", argv[i + 1]);
                    return -1;
                }
            } else {
                fprintf(stderr, "Error: --tsize option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--threads") == 0) {
            if (i + 1 < argc) {
                if (is_int(argv[i + 1])) {
                    cmd_opts[5].value = argv[i + 1];
                    i++;
                } else {
                    fprintf(stderr, "Error: --threads requires an integer value, but got '%s'.\n", argv[i + 1]);
                    return -1;
                }
            } else {
                fprintf(stderr, "Error: --threads option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--naive") == 0) {
            if (i + 1 < argc) {
                cmd_opts[6].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --naive option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--num-iter") == 0) {
            if (i + 1 < argc) {
                if (is_int(argv[i + 1])) {
                    cmd_opts[7].value = argv[i + 1];
                    i++;
                } else {
                    fprintf(stderr, "Error: --num-iter requires an integer value, but got '%s'.\n", argv[i + 1]);
                    return -1;
                }
            } else {
                fprintf(stderr, "Error: --num-iter option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cvrg-thres") == 0) {
            if (i + 1 < argc) {
                if (is_float(argv[i + 1])) {
                    cmd_opts[8].value = argv[i + 1];
                    i++;
                } else {
                    fprintf(stderr, "Error: --cvrg-thres requires an integer value, but got '%s'.\n", argv[i + 1]);
                    return -1;
                }
            } else {
                fprintf(stderr, "Error: --cvrg-thres option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--out-dir") == 0) {
            if (i + 1 < argc) {
                cmd_opts[9].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --out-dir option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--push") == 0) {
            if (i + 1 < argc) {
                cmd_opts[10].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --push option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "--dev") == 0) {
            if (i + 1 < argc) {
                cmd_opts[11].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --dev option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--relax-thres") == 0) {
            if (i + 1 < argc) {
                cmd_opts[12].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --relax option requires a value\n");
                return -1;
            }
        } else if (strcmp(argv[i], "-df") == 0 || strcmp(argv[i], "--drop-fac") == 0) {
            if (i + 1 < argc) {
                cmd_opts[13].value = argv[i + 1];
                i++;
            } else {
                fprintf(stderr, "Error: --relax option requires a value\n");
                return -1;
            }
        } else { // unknown option
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            print_help_msg(cmd_opts, n_opts);
            return -1;
        }
    }
    return 0;
}

void print_opts(opt_t *cmd_opts, int n_opts) {
    printf("Run options:\n");
    for (int i = 1; i < n_opts; i++) {
        if (cmd_opts[i].value != NULL) {
            printf("  %-15s %s\n", cmd_opts[i].flag, cmd_opts[i].value);
        } else {
            fprintf(stderr, "Error: all cmd options must be set\n");
            exit(-1);
        }
    }
}

bool check_opts(opt_t *cmd_opts, int n_opts) {
    for (int i = 1; i < n_opts; i++) {
        if (cmd_opts[i].value == NULL) {
            return false;
        }
    }
    return true;
}

double calc_rho_delt(double *prev_rho, double *curr_rho, cuint tsize) {
    double delt = 0.0;
    for (uint ti = 0; ti < tsize; ti++) {
        // delt = MAX(delt, fabs(prev_rho[ti] - curr_rho[ti]));
        delt += fabs(prev_rho[ti] - curr_rho[ti]);
    }
    return delt;
}

char* append_fn_to_dir(const char* dir, const char* fn) {
    size_t dir_len = strlen(dir);
    size_t fn_len = strlen(fn);
    if (dir_len > 0 && dir[dir_len - 1] == PATH_SEPARATOR) {
        dir_len--;
    }
    char* new_pth = (char*)malloc(dir_len + fn_len + 2);
    if (new_pth == NULL) {
        perror("Unable to allocate memory");
        return NULL;
    }
    strncpy(new_pth, dir, dir_len);
    new_pth[dir_len] = PATH_SEPARATOR;
    strncpy(new_pth + dir_len + 1, fn, fn_len + 1);
    return new_pth;
}

void drop(qinfo_t **qinfo_vec, cuint qsize, double factor) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        if (qinfo == NULL) {
            printf("qinfo[%d] is NULL\n", qi);
            exit(-1); // an error
        }

        double max_score = -1.0;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            if (tinfo->score > max_score) {
                max_score = tinfo->score;
            }
        }

        double cutoff = 1.0 / qinfo->n + (1.0 / qinfo->n * factor);
        uint ctr = 0;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            if (tinfo->rfrac < cutoff) {
                tinfo->score = 1e-2;
            } else {
                ctr++;
            }
        }

        // prevent all transcripts 
        if (ctr == 0) {
            for (uint ti = 0; ti < qinfo->n; ti++) {
                tinfo_t *tinfo = qinfo->tinfo_vec[ti];
                if (tinfo->score == max_score) {
                    tinfo->score = max_score;
                }
            }
        }

    }
}

void relax(double *rho, cuint tsize, double thres) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint ti = 0; ti < tsize; ti++) {
        if (rho[ti] < thres) {
            rho[ti] = 0.0;
        }
    }
}

void push(qinfo_t **qinfo_vec, cuint qsize) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        if (qinfo == NULL) {
            printf("qinfo[%d] is NULL\n", qi);
            exit(-1); // an error
        }
        double max_rfrac = 0.0;
        uint max_ti;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            if (tinfo->rfrac > max_rfrac) {
                max_rfrac = tinfo->rfrac;
                max_ti = ti;
            }
        }
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            tinfo->rfrac = ti == max_ti ? 1.0 : 0.0;
        }
    }
}

void push_sample(qinfo_t **qinfo_vec, cuint qsize) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint qi = 0; qi < qsize; qi++) {
        qinfo_t *qinfo = qinfo_vec[qi];
        if (qinfo == NULL) {
            printf("qinfo[%d] is NULL\n", qi);
            exit(-1); // an error
        }
        double * rfrac = malloc(qinfo->n * sizeof(double));
        double z = 0.0;
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            rfrac[ti] = tinfo->rfrac;
            z += tinfo->rfrac;
        }
        int sel_i = discrete_distn(rfrac, qinfo->n, z);
        for (uint ti = 0; ti < qinfo->n; ti++) {
            tinfo_t *tinfo = qinfo->tinfo_vec[ti];
            tinfo->rfrac = ti == sel_i ? 1.0 : 0.0;
        }
    }
}

int main(int argc, char *argv[]) {
    opt_t cmd_opts[] = {
        {"-h, --help", "help message", NULL, false}, // 0
        {"-s, --scores", "scores CSV", NULL, true}, // 1
        {"-m, --tmap", "transcripts map CSV", NULL, true}, // 2
        {"-qs, --qsize", "N of query reads", NULL, true}, // 3
        {"-ts, --tsize", "N of transcripts", NULL, true}, // 4
        {"-p, --threads", "N of threads", NULL, true}, // 5
        {"--naive", "naive mode", NULL, true}, // 6
        {"-n, --num-iter", "N of iterations", NULL, true}, // 7
        {"-c, --cvrg-thres", "convergence threshold", NULL, true}, // 8
        {"-o, --out-dir", "output directory", NULL, true}, // 9
        {"--push", "push assignments", NULL, true}, // 10
        {"--dev", "dev mode", NULL, true}, // 11
        {"-r, --relax-thres", "relax threshold", NULL, true}, // 12
        {"-df, --drop-fac", "drop factor", NULL, true} // 13
    };
    int k = 1;
    
    int n_opts = sizeof(cmd_opts) / sizeof(cmd_opts[0]);

    int parse_res = parse_args(argc, argv, cmd_opts, n_opts);

    if (parse_res == 1) {
        exit(0);
    } else if (parse_res == -1) {
        exit(1);
    }

    if (!check_opts(cmd_opts, n_opts)) {
        fprintf(stderr, "Error: all cmd options must be set\n");
        exit(-1);
    }
    cuint qsize = (cuint) atoi(cmd_opts[3].value);
    cuint tsize = (cuint) atoi(cmd_opts[4].value);
    cuint n_iter = (cuint) atoi(cmd_opts[7].value);
    double c_thres = strtod(cmd_opts[8].value, NULL);

    bool check_cvrg = true;
    if (c_thres == 0) {
        check_cvrg = false;
    }

    double r_thres = strtod(cmd_opts[12].value, NULL);
    double drop_fac = strtod(cmd_opts[13].value, NULL);

    bool is_dev = (strcmp(cmd_opts[11].value, "0") == 0) ? false : true;
    bool is_naive = (strcmp(cmd_opts[6].value, "1") == 0) ? true : false;
    bool do_push = (strcmp(cmd_opts[10].value, "1") == 0) ? true : false;
    bool do_drop = (drop_fac == 0.0) ? false : true;
    bool do_relax = (r_thres == 0.0) ? false : true;

    cuint n_th = (cuint) atoi(cmd_opts[5].value);

    #ifdef _OPENMP
    omp_set_num_threads(n_th);
    #endif

    qinfo_t ** qinfo_vec = load_qt_info(cmd_opts[1].value, qsize, 200);
    char ** tnames = (char**) malloc(tsize * sizeof(char*));
    hmap_t * tmap = load_tmap(cmd_opts[2].value, tnames);

    FILE *train_file;
    FILE *delt_file;
    if (is_dev) {
        char * fn = append_fn_to_dir(cmd_opts[9].value, "train.csv");
        train_file = fopen(fn, "w");
        fprintf(train_file, "i");
        for (uint ti = 0; ti < tsize; ti++) {
            char * tname = tnames[ti];
            fprintf(train_file, ",%s", tname);
        }
        fprintf(train_file, "\n");
        char * fn_2 = append_fn_to_dir(cmd_opts[9].value, "delt.txt");
        delt_file = fopen(fn_2, "w");
        fprintf(delt_file, "i,delt\n");
    }
    double *rho = init(qsize, tsize);
    double *rho_prev = (double *) malloc(tsize * sizeof(double));
    double rho_delt;
    uint ctr = 0;
    for (uint i = 0; i < n_iter; i++) {
        step_e(qinfo_vec, rho, tmap, qsize, is_naive);

        if (i == 0 && do_drop) {
            drop(qinfo_vec, qsize, drop_fac);
            step_e(qinfo_vec, rho, tmap, qsize, is_naive);
        }

        memcpy(rho_prev, rho, tsize * sizeof(double));
        memset(rho, 0, tsize * sizeof(double));
        step_m(qinfo_vec, rho, tmap, qsize);

        if (r_thres > 0.0) {
            relax(rho, tsize, r_thres);
        }
        
        rho_delt = calc_rho_delt(rho_prev, rho, tsize);

        if (check_cvrg && rho_delt < c_thres) {
            break;
        }
        
        if (is_dev) {
            fprintf(train_file, "%d", i);
            for (uint ti = 0; ti < tsize; ti++) {
                fprintf(train_file, ",%f", rho[ti]);
            }
            fprintf(train_file, "\n");
            fprintf(delt_file, "%d,%f\n", i, rho_delt);
        }
        ctr++;
    }
    
    printf("N of EM iterations: %d\n", ctr);
    printf("Final delta: %f\n", rho_delt);

    if (is_dev) {
        fclose(train_file);
        fclose(delt_file);
    }

    if (do_push) {
        srand(time(NULL));
        push_sample(qinfo_vec, qsize);
        // push(qinfo_vec, qsize);
        memcpy(rho_prev, rho, tsize * sizeof(double));
        memset(rho, 0, tsize * sizeof(double));
        step_m(qinfo_vec, rho, tmap, qsize);
        rho_delt = calc_rho_delt(rho_prev, rho, tsize);
        printf("Delta following push: %f\n", rho_delt);
    }

    char * asgn_fn = append_fn_to_dir(cmd_opts[9].value, "assignments.csv");
    char * unasgn_fn = append_fn_to_dir(cmd_opts[9].value, "unassigned.txt");
    char * abnd_fn = append_fn_to_dir(cmd_opts[9].value, "abundances.csv");

    write_results((const qinfo_t **) qinfo_vec, qsize, tsize, \
                rho, asgn_fn, unasgn_fn, abnd_fn, tnames);

    hmap_free(tmap);
    free_qt_info(qinfo_vec, qsize, tsize);
}