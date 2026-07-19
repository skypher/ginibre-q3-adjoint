#ifdef USE_GMP
#include <gmpxx.h>
#else
#include <boost/multiprecision/cpp_int.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef USE_GMP
using BigInt = mpz_class;
#else
using BigInt = boost::multiprecision::cpp_int;
#endif

using Key = std::vector<int>;

struct KeyHash {
    std::size_t operator()(const Key& v) const {
        std::uint64_t h = 1469598103934665603ull;
        for (int x : v) {
            std::uint64_t y = static_cast<std::uint64_t>(static_cast<std::int64_t>(x));
            h ^= y + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

struct Weight {
    Key w;
    int mult;
};

struct ReflectResult {
    bool wall = false;
    int sign = 1;
    Key dom;
};

std::vector<std::vector<int>> cartan_B(int n) {
    std::vector<std::vector<int>> A(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) A[i][i] = 2;
    for (int i = 0; i < n - 1; ++i) {
        if (i < n - 2) {
            A[i][i + 1] = -1;
            A[i + 1][i] = -1;
        } else {
            A[i][i + 1] = -2;
            A[i + 1][i] = -1;
        }
    }
    return A;
}

std::vector<std::vector<int>> cartan_C(int n) {
    std::vector<std::vector<int>> A(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) A[i][i] = 2;
    for (int i = 0; i < n - 1; ++i) {
        if (i < n - 2) {
            A[i][i + 1] = -1;
            A[i + 1][i] = -1;
        } else {
            A[i][i + 1] = -1;
            A[i + 1][i] = -2;
        }
    }
    return A;
}

std::vector<std::vector<int>> cartan_D(int n) {
    if (n < 4) throw std::runtime_error("D rank must be at least 4");
    std::vector<std::vector<int>> A(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) A[i][i] = 2;
    A[0][2] = -1;
    A[2][0] = -1;
    A[1][2] = -1;
    A[2][1] = -1;
    for (int i = 2; i < n - 1; ++i) {
        A[i][i + 1] = -1;
        A[i + 1][i] = -1;
    }
    return A;
}

std::vector<std::vector<int>> cartan_matrix(char family, int rank) {
    if (family == 'B') return cartan_B(rank);
    if (family == 'C') return cartan_C(rank);
    if (family == 'D') return cartan_D(rank);
    throw std::runtime_error("unknown family");
}

int height(const Key& v) {
    int s = 0;
    for (int x : v) s += x;
    return s;
}

std::vector<Key> positive_roots_simple(const std::vector<std::vector<int>>& cartan) {
    int r = static_cast<int>(cartan.size());
    std::unordered_set<Key, KeyHash> roots;
    std::vector<Key> frontier;
    for (int i = 0; i < r; ++i) {
        Key e(r, 0);
        e[i] = 1;
        roots.insert(e);
        frontier.push_back(e);
    }
    while (!frontier.empty()) {
        std::vector<Key> next;
        for (const Key& beta : frontier) {
            for (int i = 0; i < r; ++i) {
                int pairing = 0;
                for (int j = 0; j < r; ++j) pairing += beta[j] * cartan[j][i];
                int p = 0;
                while (true) {
                    Key cand = beta;
                    cand[i] -= (p + 1);
                    if (roots.find(cand) != roots.end()) {
                        ++p;
                    } else {
                        break;
                    }
                }
                int q_max = p - pairing;
                if (q_max > 0) {
                    Key cand = beta;
                    cand[i] += 1;
                    if (roots.insert(cand).second) next.push_back(std::move(cand));
                }
            }
        }
        frontier = std::move(next);
    }
    std::vector<Key> out(roots.begin(), roots.end());
    std::sort(out.begin(), out.end(), [](const Key& a, const Key& b) {
        int ha = height(a), hb = height(b);
        if (ha != hb) return ha < hb;
        return a < b;
    });
    return out;
}

Key simple_to_fundamental(const Key& beta, const std::vector<std::vector<int>>& cartan) {
    int r = static_cast<int>(cartan.size());
    Key fw(r, 0);
    for (int k = 0; k < r; ++k) {
        int s = 0;
        for (int i = 0; i < r; ++i) s += beta[i] * cartan[i][k];
        fw[k] = s;
    }
    return fw;
}

std::vector<Weight> adjoint_weights(const std::vector<std::vector<int>>& cartan) {
    int r = static_cast<int>(cartan.size());
    std::vector<Weight> weights;
    for (const Key& beta : positive_roots_simple(cartan)) {
        Key fw = simple_to_fundamental(beta, cartan);
        weights.push_back({fw, 1});
        for (int& x : fw) x = -x;
        weights.push_back({fw, 1});
    }
    weights.push_back({Key(r, 0), r});
    return weights;
}

class LieGroup {
   public:
    explicit LieGroup(std::vector<std::vector<int>> cartan, std::size_t reflect_cache_limit)
        : cartan_(std::move(cartan)),
          rank_(static_cast<int>(cartan_.size())),
          adj_(adjoint_weights(cartan_)),
          reflect_cache_limit_(reflect_cache_limit) {}

    ReflectResult dominant_reflect(const Key& mu) {
        if (reflect_cache_limit_ != 0) {
            auto it = cache_.find(mu);
            if (it != cache_.end()) return it->second;
        }

        Key n = mu;
        for (int& x : n) ++x;
        int sign = 1;
        for (int step = 0; step < 10000; ++step) {
            int bad = -1;
            for (int i = 0; i < rank_; ++i) {
                if (n[i] < 0) {
                    bad = i;
                    break;
                }
                if (n[i] == 0) {
                    ReflectResult res;
                    res.wall = true;
                    maybe_cache(mu, res);
                    return res;
                }
            }
            if (bad < 0) break;
            int ni = n[bad];
            for (int k = 0; k < rank_; ++k) n[k] -= ni * cartan_[bad][k];
            sign = -sign;
        }
        for (int& x : n) --x;
        ReflectResult res;
        res.sign = sign;
        res.dom = std::move(n);
        maybe_cache(mu, res);
        return res;
    }

    std::unordered_map<Key, BigInt, KeyHash> tensor_with_adj(const std::unordered_map<Key, BigInt, KeyHash>& V) {
        std::unordered_map<Key, BigInt, KeyHash> W;
        W.reserve(V.size() * 2);
        Key target(rank_);
        for (const auto& item : V) {
            const Key& mu = item.first;
            const BigInt& coeff = item.second;
            if (coeff == 0) continue;
            for (const Weight& eta : adj_) {
                for (int i = 0; i < rank_; ++i) target[i] = mu[i] + eta.w[i];
                ReflectResult res = dominant_reflect(target);
                if (res.wall) continue;
                if (res.sign > 0) {
                    W[res.dom] += coeff * eta.mult;
                } else {
                    W[res.dom] -= coeff * eta.mult;
                }
            }
        }
        for (auto it = W.begin(); it != W.end();) {
            if (it->second == 0) {
                it = W.erase(it);
            } else {
                ++it;
            }
        }
        return W;
    }

    std::vector<BigInt> compute_moments(int max_k, bool verbose, bool clear_cache_each_moment) {
        Key zero(rank_, 0);
        std::unordered_map<Key, BigInt, KeyHash> V;
        V.emplace(zero, BigInt(1));
        std::vector<BigInt> moments;
        moments.push_back(BigInt(1));
        for (int k = 1; k <= max_k; ++k) {
            V = tensor_with_adj(V);
            auto it = V.find(zero);
            moments.push_back(it == V.end() ? BigInt(0) : it->second);
            if (verbose) {
                std::cout << "  m_" << k << " = " << moments.back()
                          << "  (|support| = " << V.size()
                          << ", |cache| = " << cache_.size() << ")" << std::endl;
            }
            if (clear_cache_each_moment && k < max_k) {
                std::unordered_map<Key, ReflectResult, KeyHash>().swap(cache_);
#ifdef __GLIBC__
                malloc_trim(0);
#endif
            }
        }
        return moments;
    }

   private:
    void maybe_cache(const Key& mu, const ReflectResult& res) {
        if (cache_.size() < reflect_cache_limit_) cache_.emplace(mu, res);
    }

    std::vector<std::vector<int>> cartan_;
    int rank_;
    std::vector<Weight> adj_;
    std::unordered_map<Key, ReflectResult, KeyHash> cache_;
    std::size_t reflect_cache_limit_;
};

BigInt binom_int(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    BigInt ans = 1;
    for (int i = 1; i <= k; ++i) {
        ans *= (n - k + i);
        ans /= i;
    }
    return ans;
}

BigInt q3(const std::vector<BigInt>& moments, int n) {
    BigInt total = 0;
    for (int k = 0; k <= n; ++k) {
        total += binom_int(n, k) * (moments[k + 2] * moments[n - k] - moments[k + 1] * moments[n - k + 1]);
    }
    return 2 * total;
}

BigInt chain_diff(const std::vector<BigInt>& moments, int m) {
    return q3(moments, 2 * m + 3) - 4 * q3(moments, 2 * m + 1);
}

struct Args {
    int max_chain = 4;
    std::string families = "B,C,D";
    int rank_min = -1;
    int rank_cap = -1;
    bool moment_progress = false;
    bool quick_exit = false;
    bool clear_cache_each_moment = false;
    std::size_t reflect_cache_limit = std::numeric_limits<std::size_t>::max();
};

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + flag);
            return argv[++i];
        };
        if (a == "--max-chain") args.max_chain = std::stoi(need_value(a));
        else if (a == "--families") args.families = need_value(a);
        else if (a == "--rank-min") args.rank_min = std::stoi(need_value(a));
        else if (a == "--rank-cap") args.rank_cap = std::stoi(need_value(a));
        else if (a == "--moment-progress") args.moment_progress = true;
        else if (a == "--quick-exit") args.quick_exit = true;
        else if (a == "--clear-cache-each-moment") args.clear_cache_each_moment = true;
        else if (a == "--reflect-cache-limit") args.reflect_cache_limit = static_cast<std::size_t>(std::stoull(need_value(a)));
        else throw std::runtime_error("unknown argument " + a);
    }
    return args;
}

std::string cache_limit_label(std::size_t limit) {
    if (limit == std::numeric_limits<std::size_t>::max()) return "unlimited";
    return std::to_string(limit);
}

std::vector<char> parse_families(const std::string& s) {
    std::vector<char> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(tok[0]);
    }
    return out;
}

std::vector<int> nonstable_ranks(char family, int max_chain) {
    int lo = family == 'D' ? 4 : 2;
    int hi = 0;
    if (family == 'B' || family == 'C') hi = max_chain + 1;
    else if (family == 'D') hi = 2 * max_chain + 5;
    else throw std::runtime_error("unknown family");
    std::vector<int> ranks;
    for (int r = lo; r <= hi; ++r) ranks.push_back(r);
    return ranks;
}

int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);
        int moment_max = 2 * args.max_chain + 5;
        std::cout << "Exact classical Chain bridge C++: max_chain=" << args.max_chain
                  << ", moment_max=" << moment_max
                  << ", families=" << args.families
                  << ", rank_min=" << args.rank_min
                  << ", rank_cap=" << args.rank_cap
                  << ", reflect_cache_limit=" << cache_limit_label(args.reflect_cache_limit) << std::endl;

        bool all_ok = true;
        struct Row {
            std::string group;
            int rank;
            BigInt min_diff;
            BigInt last_diff;
            BigInt last_q3;
            double seconds;
        };
        std::vector<Row> rows;

        for (char family : parse_families(args.families)) {
            for (int rank : nonstable_ranks(family, args.max_chain)) {
                if (args.rank_min >= 0 && rank < args.rank_min) continue;
                if (args.rank_cap >= 0 && rank > args.rank_cap) continue;
                std::string group = std::string(1, family) + std::to_string(rank);
                auto t0 = std::chrono::steady_clock::now();
                LieGroup G(cartan_matrix(family, rank), args.reflect_cache_limit);
                std::vector<BigInt> moments =
                    G.compute_moments(moment_max, args.moment_progress, args.clear_cache_each_moment);
                std::vector<BigInt> diffs;
                std::vector<BigInt> qs;
                for (int m = 0; m <= args.max_chain; ++m) diffs.push_back(chain_diff(moments, m));
                for (int m = 0; m <= args.max_chain + 1; ++m) qs.push_back(q3(moments, 2 * m + 1));
                BigInt min_diff = diffs[0];
                for (const BigInt& d : diffs) {
                    if (d < min_diff) min_diff = d;
                }
                bool ok = true;
                for (const BigInt& d : diffs) ok = ok && d >= 0;
                for (const BigInt& q : qs) ok = ok && q >= 0;
                all_ok = all_ok && ok;
                auto t1 = std::chrono::steady_clock::now();
                double seconds = std::chrono::duration<double>(t1 - t0).count();
                rows.push_back({group, rank, min_diff, diffs.back(), qs.back(), seconds});
                std::cout << group << ": min_diff_m0.." << args.max_chain << "=" << min_diff
                          << " diff_m" << args.max_chain << "=" << diffs.back()
                          << " Q3(" << (2 * args.max_chain + 3) << ")=" << qs.back()
                          << " m0..6=[";
                for (int i = 0; i <= 6 && i < static_cast<int>(moments.size()); ++i) {
                    if (i) std::cout << ", ";
                    std::cout << moments[i];
                }
                std::cout << "] dt=" << seconds << "s " << (ok ? "OK" : "FAIL") << std::endl;
                if (args.quick_exit && args.rank_min == rank && args.rank_cap == rank) {
                    std::cout << "\nSummary table:\n";
                    std::cout << "group,rank,min_chain_diff,last_chain_diff,last_odd_Q3,seconds\n";
                    const Row& row = rows.back();
                    std::cout << row.group << "," << row.rank << "," << row.min_diff << ","
                              << row.last_diff << "," << row.last_q3 << "," << row.seconds << "\n";
                    std::cout << "\nRESULT: " << (all_ok ? "ALL PASS" : "FAIL") << std::endl;
                    std::cout.flush();
                    std::cerr.flush();
                    std::quick_exit(all_ok ? 0 : 1);
                }
            }
        }

        std::cout << "\nSummary table:\n";
        std::cout << "group,rank,min_chain_diff,last_chain_diff,last_odd_Q3,seconds\n";
        for (const Row& row : rows) {
            std::cout << row.group << "," << row.rank << "," << row.min_diff << ","
                      << row.last_diff << "," << row.last_q3 << "," << row.seconds << "\n";
        }
        std::cout << "\nRESULT: " << (all_ok ? "ALL PASS" : "FAIL") << std::endl;
        if (args.quick_exit) {
            std::cout.flush();
            std::cerr.flush();
            std::quick_exit(all_ok ? 0 : 1);
        }
        return all_ok ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << std::endl;
        return 2;
    }
}
