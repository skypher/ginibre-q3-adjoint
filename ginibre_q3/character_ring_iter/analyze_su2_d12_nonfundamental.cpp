#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Minus = std::array<int, 2>;
using Plus = std::array<int, 5>;

template<class F>
void outputs(int k, int a, int b, F f) {
    const int upper = std::min(a+b, 2*k-a-b);
    for (int c = std::abs(a-b); c <= upper; c += 2) f(c);
}

std::int64_t mult(int k, const std::vector<int>& labels) {
    std::vector<std::int64_t> current(static_cast<std::size_t>(k+1));
    std::vector<std::int64_t> next(static_cast<std::size_t>(k+1));
    current[0] = 1;
    for (int p : labels) {
        std::fill(next.begin(), next.end(), 0);
        for (int x = 0; x <= k; ++x) {
            const auto xi = static_cast<std::size_t>(x);
            if (current[xi] == 0) continue;
            outputs(k, x, p, [&](int y) {
                next[static_cast<std::size_t>(y)] += current[xi];
            });
        }
        current.swap(next);
    }
    return current[0];
}

bool disjoint(const Minus& minus, const Plus& plus) {
    for (int q : minus) for (int p : plus) if (q == p) return false;
    return true;
}

struct Stats {
    std::int64_t n{}, peq{}, u1{}, u2{}, t{}, d{};
    int active{};
};

Stats statistics(int k, const Minus& minus, const Plus& plus) {
    Stats s;
    s.n = mult(k, {
        minus[0], minus[1], plus[0], plus[1], plus[2], plus[3], plus[4]
    });
    if (minus[0] == minus[1]) {
        s.peq += mult(k, {plus[0], plus[1], plus[2], plus[3], plus[4]});
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        std::vector<int> rest;
        for (std::size_t j = 0; j < plus.size(); ++j) {
            if (j != i) rest.push_back(plus[j]);
        }
        s.u1 += mult(k, {minus[0], minus[1], plus[i]}) * mult(k, rest);
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
      for (std::size_t j = i+1U; j < plus.size(); ++j) {
        std::vector<int> complement;
        std::vector<int> positive_rest{minus[0], minus[1]};
        for (std::size_t h = 0; h < plus.size(); ++h) {
            if (h == i || h == j) continue;
            complement.push_back(plus[h]);
            positive_rest.push_back(plus[h]);
        }
        if (plus[i] == plus[j]) s.peq += mult(k, positive_rest);
        s.u2 += mult(k, {minus[0],minus[1],plus[i],plus[j]})
              * mult(k, complement);
        for (std::size_t o = 0; o < minus.size(); ++o) {
            const auto triple = mult(k, {minus[o],plus[i],plus[j]});
            if (triple == 0) continue;
            auto block = complement;
            block.push_back(minus[1U-o]);
            const auto cut = triple * mult(k, block);
            if (cut == 0) continue;
            ++s.active;
            s.t += cut;
            s.d = std::max(s.d, cut);
        }
      }
    }
    return s;
}

struct Case {
    int k{};
    Minus minus{};
    Plus plus{};
    Stats s{};
};

std::string signature(const Case& z) {
    std::string out = "k=" + std::to_string(z.k) + " m=";
    for (int x : z.minus) out += std::to_string(x) + ",";
    out += " p=";
    for (int x : z.plus) out += std::to_string(x) + ",";
    return out;
}

void print_case(const char* tag, const Case& z) {
    const auto raw = z.s.n + z.s.peq - z.s.t;
    std::cout << tag << " k=" << z.k << " minus=["
              << z.minus[0] << ',' << z.minus[1] << "] plus=[";
    for (std::size_t i = 0; i < z.plus.size(); ++i) {
        std::cout << (i ? "," : "") << z.plus[i];
    }
    std::cout << "] d=" << z.s.d << " c=" << z.s.active
              << " N=" << z.s.n << " P=" << z.s.peq
              << " U1=" << z.s.u1 << " U2=" << z.s.u2
              << " T=" << z.s.t << " raw=" << raw
              << " paid=" << raw+z.s.u1+z.s.u2 << '\n';
}

std::string cap3_name(int k, int x) {
    if (x == 1) return "1";
    if (x == 2) return "2";
    if (x == k-2) return "h";
    if (x == k-1) return "e";
    if (x == k) return "J";
    return "?";
}

void print_cap3_case(const Case& z) {
    const auto raw = z.s.n + z.s.peq - z.s.t;
    std::cout << "m=" << cap3_name(z.k,z.minus[0]) << ','
              << cap3_name(z.k,z.minus[1]) << " p=";
    for (std::size_t i = 0; i < z.plus.size(); ++i) {
        std::cout << (i ? "," : "") << cap3_name(z.k,z.plus[i]);
    }
    std::cout << " d=" << z.s.d << " N=" << z.s.n
              << " P=" << z.s.peq << " T=" << z.s.t
              << " raw=" << raw << '\n';
}

void print_detail(int k, const Minus& minus, const Plus& plus) {
    for (std::size_t i = 0; i < plus.size(); ++i) {
      for (std::size_t j = i+1U; j < plus.size(); ++j) {
        std::vector<int> complement;
        for (std::size_t h = 0; h < plus.size(); ++h) {
            if (h != i && h != j) complement.push_back(plus[h]);
        }
        for (std::size_t o = 0; o < minus.size(); ++o) {
            const auto triple = mult(k, {minus[o],plus[i],plus[j]});
            auto block = complement;
            block.push_back(minus[1U-o]);
            const auto rank = mult(k, block);
            if (triple*rank == 0) continue;
            std::cout << "cut o=" << o << " pair=(" << i << ',' << j
                      << ") labels=(" << plus[i] << ',' << plus[j]
                      << ") triple=" << triple << " rank=" << rank
                      << " value=" << triple*rank << '\n';
        }
      }
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
        std::vector<int> rest;
        for (std::size_t h = 0; h < plus.size(); ++h) {
            if (h != i) rest.push_back(plus[h]);
        }
        const auto left=mult(k,{minus[0],minus[1],plus[i]});
        const auto right=mult(k,rest);
        if (left*right != 0) {
            std::cout << "U1 i=" << i << " label=" << plus[i]
                      << " factors=(" << left << ',' << right
                      << ") value=" << left*right << '\n';
        }
    }
    for (std::size_t i = 0; i < plus.size(); ++i) {
      for (std::size_t j = i+1U; j < plus.size(); ++j) {
        std::vector<int> rest;
        for (std::size_t h = 0; h < plus.size(); ++h) {
            if (h != i && h != j) rest.push_back(plus[h]);
        }
        const auto left=mult(k,{minus[0],minus[1],plus[i],plus[j]});
        const auto right=mult(k,rest);
        if (left*right != 0) {
            std::cout << "U2 pair=(" << i << ',' << j << ") labels=("
                      << plus[i] << ',' << plus[j] << ") factors=("
                      << left << ',' << right << ") value="
                      << left*right << '\n';
        }
      }
    }
}

void record_case(
    int k, Minus minus, Plus plus,
    std::vector<Case>& local_deficits,
    std::array<Case,3>& local_paid_min,
    std::array<bool,3>& local_paid_set
) {
    if (minus == Minus{1,1} || !disjoint(minus, plus)) return;
    const auto s = statistics(k, minus, plus);
    if (s.d < 1 || s.d > 2) return;
    const auto index = static_cast<std::size_t>(s.d);
    const auto raw = s.n+s.peq-s.t;
    const auto paid = raw+s.u1+s.u2;
    Case z{k,minus,plus,s};
    if (raw < 0) local_deficits.push_back(z);
    if (!local_paid_set[index]
        || paid < local_paid_min[index].s.n
                 +local_paid_min[index].s.peq
                 +local_paid_min[index].s.u1
                 +local_paid_min[index].s.u2
                 -local_paid_min[index].s.t) {
        local_paid_set[index] = true;
        local_paid_min[index] = z;
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 3
        && (std::string(argv[1]) == "--cap3-table"
            || std::string(argv[1]) == "--cap3-full"
            || std::string(argv[1]) == "--cap3-summary")) {
        const int k = std::atoi(argv[2]);
        if (k < 6) return EXIT_FAILURE;
        const bool full = std::string(argv[1]) == "--cap3-full";
        const bool summary = std::string(argv[1]) == "--cap3-summary";
        const std::array<int,5> labels{1,2,k-2,k-1,k};
        for (std::size_t qi=0; qi<labels.size(); ++qi) {
          for (std::size_t ai=qi; ai<labels.size(); ++ai) {
            const Minus minus{labels[qi],labels[ai]};
            if (minus==Minus{1,1}) continue;
            bool set=false;
            Case best;
            int total=0, zero=0, small=0, large=0, deficits=0;
            for (std::size_t i0=0;i0<labels.size();++i0)
            for (std::size_t i1=i0;i1<labels.size();++i1)
            for (std::size_t i2=i1;i2<labels.size();++i2)
            for (std::size_t i3=i2;i3<labels.size();++i3)
            for (std::size_t i4=i3;i4<labels.size();++i4) {
                const Plus plus{
                    labels[i0],labels[i1],labels[i2],labels[i3],labels[i4]
                };
                if (!disjoint(minus,plus)) continue;
                ++total;
                const auto s=statistics(k,minus,plus);
                if (s.d == 0) {
                    ++zero;
                    continue;
                }
                if (s.d > 2) {
                    ++large;
                    continue;
                }
                ++small;
                if (full) print_cap3_case(Case{k,minus,plus,s});
                const auto raw=s.n+s.peq-s.t;
                if (raw < 0) ++deficits;
                if (!set || raw<best.s.n+best.s.peq-best.s.t) {
                    set=true;
                    best=Case{k,minus,plus,s};
                }
            }
            if (summary) {
                std::cout << "summary m=" << cap3_name(k,minus[0]) << ','
                          << cap3_name(k,minus[1])
                          << " total=" << total << " zero=" << zero
                          << " small=" << small << " large=" << large
                          << " deficits=" << deficits;
                if (set) {
                    std::cout << " min=" << best.s.n+best.s.peq-best.s.t
                              << " witness=";
                    for (std::size_t i=0;i<best.plus.size();++i) {
                        std::cout << (i ? "," : "")
                                  << cap3_name(k,best.plus[i]);
                    }
                }
                std::cout << '\n';
            } else if (!full && set) {
                print_case("cap3-min",best);
            }
          }
        }
        return EXIT_SUCCESS;
    }
    if (argc == 10 && std::string(argv[1]) == "--case") {
        const int k=std::atoi(argv[2]);
        Minus minus{std::atoi(argv[3]),std::atoi(argv[4])};
        Plus plus{
            std::atoi(argv[5]),std::atoi(argv[6]),std::atoi(argv[7]),
            std::atoi(argv[8]),std::atoi(argv[9])
        };
        std::sort(minus.begin(),minus.end());
        std::sort(plus.begin(),plus.end());
        print_case("case",Case{k,minus,plus,statistics(k,minus,plus)});
        print_detail(k,minus,plus);
        return EXIT_SUCCESS;
    }
    if (argc != 4) {
        std::cerr << "usage: analyze_su2_d12_nonfundamental"
                     " --random MAX_K SAMPLES\n"
                     "   or: analyze_su2_d12_nonfundamental"
                     " --exact MIN_K MAX_K\n"
                     "   or: analyze_su2_d12_nonfundamental"
                     " --case K Q A P1 P2 P3 P4 P5\n";
        return EXIT_FAILURE;
    }
    const std::string mode = argv[1];
    const int first = std::atoi(argv[2]);
    const std::uint64_t last = std::strtoull(argv[3], nullptr, 10);
    std::vector<Case> deficits;
    std::array<Case,3> paid_min{};
    std::array<bool,3> paid_set{};
    std::mutex merge_mutex;

    if (mode == "--random") {
        const int max_k = first;
        const std::uint64_t samples = last;
#pragma omp parallel
        {
            std::vector<Case> local_deficits;
            std::array<Case,3> local_paid_min{};
            std::array<bool,3> local_paid_set{};
#pragma omp for schedule(static)
            for (std::uint64_t sample = 0; sample < samples; ++sample) {
                std::uint64_t state = sample + 0x9e3779b97f4a7c15ULL;
                auto random = [&state]() {
                    state += 0x9e3779b97f4a7c15ULL;
                    std::uint64_t z = state;
                    z = (z^(z>>30U))*0xbf58476d1ce4e5b9ULL;
                    z = (z^(z>>27U))*0x94d049bb133111ebULL;
                    return z^(z>>31U);
                };
                const int k = 4 + static_cast<int>(
                    random()%static_cast<std::uint64_t>(max_k-3)
                );
                Minus minus{};
                Plus plus{};
                for (int& x : minus) {
                    x = 1+static_cast<int>(
                        random()%static_cast<std::uint64_t>(k)
                    );
                }
                for (int& x : plus) {
                    x = 1+static_cast<int>(
                        random()%static_cast<std::uint64_t>(k)
                    );
                }
                std::sort(minus.begin(),minus.end());
                std::sort(plus.begin(),plus.end());
                record_case(
                    k,minus,plus,local_deficits,local_paid_min,local_paid_set
                );
            }
            std::lock_guard<std::mutex> lock(merge_mutex);
            deficits.insert(
                deficits.end(),local_deficits.begin(),local_deficits.end()
            );
            for (std::size_t d = 1; d <= 2; ++d) {
                if (!local_paid_set[d]) continue;
                std::vector<Case> sink;
                record_case(
                    local_paid_min[d].k,local_paid_min[d].minus,
                    local_paid_min[d].plus,sink,paid_min,paid_set
                );
            }
        }
    } else if (mode == "--exact") {
        const int min_k = first;
        const int max_k = static_cast<int>(last);
        for (int k = min_k; k <= max_k; ++k) {
#pragma omp parallel
            {
                std::vector<Case> local_deficits;
                std::array<Case,3> local_paid_min{};
                std::array<bool,3> local_paid_set{};
#pragma omp for schedule(dynamic)
                for (int q = 1; q <= k; ++q) {
                  for (int a = q; a <= k; ++a) {
                    const Minus minus{q,a};
                    for (int p1=1;p1<=k;++p1)
                    for (int p2=p1;p2<=k;++p2)
                    for (int p3=p2;p3<=k;++p3)
                    for (int p4=p3;p4<=k;++p4)
                    for (int p5=p4;p5<=k;++p5) {
                        const Plus plus{p1,p2,p3,p4,p5};
                        record_case(
                            k,minus,plus,local_deficits,
                            local_paid_min,local_paid_set
                        );
                    }
                  }
                }
                std::lock_guard<std::mutex> lock(merge_mutex);
                deficits.insert(
                    deficits.end(),local_deficits.begin(),local_deficits.end()
                );
                for (std::size_t d = 1; d <= 2; ++d) {
                    if (!local_paid_set[d]) continue;
                    std::vector<Case> sink;
                    record_case(
                        local_paid_min[d].k,local_paid_min[d].minus,
                        local_paid_min[d].plus,sink,paid_min,paid_set
                    );
                }
            }
            std::cerr << "finished k=" << k
                      << " deficits=" << deficits.size() << '\n';
        }
    } else {
        std::cerr << "unknown mode\n";
        return EXIT_FAILURE;
    }

    std::sort(deficits.begin(),deficits.end(),[](const Case& x,const Case& y){
        const auto xm=x.s.n+x.s.peq-x.s.t;
        const auto ym=y.s.n+y.s.peq-y.s.t;
        return std::tie(xm,x.k,x.minus,x.plus)<std::tie(ym,y.k,y.minus,y.plus);
    });
    std::set<std::string> seen;
    std::size_t printed=0;
    for (const auto& z: deficits) {
        if (!seen.insert(signature(z)).second) continue;
        if (printed++<1000U) print_case("deficit",z);
    }
    for (std::size_t d=1;d<=2;++d) {
        if (paid_set[d]) print_case("paid-min",paid_min[d]);
    }
    std::cout << "deficit-count=" << deficits.size()
              << " unique-printed=" << std::min<std::size_t>(printed,1000U)
              << '\n';
}
