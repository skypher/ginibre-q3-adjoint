#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using Integer = boost::multiprecision::cpp_int;

Integer at(const std::vector<Integer>& profile, int index) {
    if (index < 0 || index >= static_cast<int>(profile.size())) {
        return 0;
    }
    return profile[static_cast<std::size_t>(index)];
}

std::vector<Integer> parse_profile(const std::string& text) {
    std::vector<Integer> result;
    std::stringstream input{text};
    std::string token;
    while (std::getline(input, token, ',')) {
        if (token.empty()) {
            throw std::invalid_argument("empty profile entry");
        }
        std::stringstream entry{token};
        Integer value;
        entry >> value;
        if (!entry || !entry.eof() || value < 0) {
            throw std::invalid_argument(
                "profile entries must be nonnegative integers"
            );
        }
        result.push_back(value);
    }
    if (result.empty()) {
        throw std::invalid_argument("profile must be nonempty");
    }
    return result;
}

std::vector<Integer> transform(
    const std::vector<Integer>& profile,
    int label
) {
    std::vector<Integer> result(
        profile.size() + static_cast<std::size_t>(label)
    );
    for (int target = 0; target < static_cast<int>(result.size());
         ++target) {
        for (int source = std::abs(target - label);
             source <= target + label;
             ++source) {
            result[static_cast<std::size_t>(target)] +=
                at(profile, source);
        }
    }
    return result;
}

Integer inner(
    const std::vector<Integer>& left,
    const std::vector<Integer>& right,
    int cutoff
) {
    Integer result = 0;
    const int size = std::max(
        static_cast<int>(left.size()),
        static_cast<int>(right.size())
    );
    for (int index = cutoff; index < size; ++index) {
        result += at(left, index) * at(right, index);
    }
    return result;
}

bool log_concave(const std::vector<Integer>& profile) {
    bool seen_positive = false;
    bool ended = false;
    for (std::size_t index = 0; index < profile.size(); ++index) {
        if (profile[index] == 0) {
            if (seen_positive) {
                ended = true;
            }
        } else {
            if (ended) {
                return false;
            }
            seen_positive = true;
        }
        if (
            index > 0 && index + 1 < profile.size()
            && profile[index] * profile[index]
                < profile[index - 1] * profile[index + 1]
        ) {
            return false;
        }
    }
    return seen_positive;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr
            << "usage: analyze_su2_log_concave_tail_profile"
            << " p0,p1,... max_label\n";
        return EXIT_FAILURE;
    }
    try {
        const std::vector<Integer> profile = parse_profile(argv[1]);
        const int max_label = std::stoi(argv[2]);
        if (max_label < 1 || !log_concave(profile)) {
            throw std::invalid_argument(
                "profile must be log-concave and max_label positive"
            );
        }
        std::vector<std::vector<Integer>> images;
        images.push_back(profile);
        for (int label = 1; label <= max_label; ++label) {
            images.push_back(transform(profile, label));
        }
        Integer minimum = 0;
        bool initialized = false;
        int minimum_q = 0;
        int minimum_a = 0;
        int minimum_cutoff = 0;
        for (int q = 1; q <= max_label; ++q) {
            for (int a = 1; a <= max_label; ++a) {
                const int size = std::max(
                    static_cast<int>(
                        images[static_cast<std::size_t>(q)].size()
                    ),
                    static_cast<int>(
                        images[static_cast<std::size_t>(a)].size()
                    )
                );
                for (int cutoff = 0; cutoff <= size; ++cutoff) {
                    const Integer determinant =
                        inner(profile, profile, cutoff)
                            * inner(
                                images[static_cast<std::size_t>(q)],
                                images[static_cast<std::size_t>(a)],
                                cutoff
                            )
                        - inner(
                              profile,
                              images[static_cast<std::size_t>(q)],
                              cutoff
                          )
                            * inner(
                                profile,
                                images[static_cast<std::size_t>(a)],
                                cutoff
                            );
                    if (!initialized || determinant < minimum) {
                        initialized = true;
                        minimum = determinant;
                        minimum_q = q;
                        minimum_a = a;
                        minimum_cutoff = cutoff;
                    }
                }
            }
        }
        Integer full_minimum = 0;
        bool full_initialized = false;
        int full_r = 0;
        int full_s = 0;
        int full_c = 0;
        int full_d = 0;
        int full_cutoff = 0;
        for (int r = 0; r <= max_label; ++r) {
            for (int s = r + 1; s <= max_label; ++s) {
                for (int c = 0; c <= max_label; ++c) {
                    for (int d = c + 1; d <= max_label; ++d) {
                        const int size = std::max(
                            {
                                static_cast<int>(
                                    images[
                                        static_cast<std::size_t>(r)
                                    ].size()
                                ),
                                static_cast<int>(
                                    images[
                                        static_cast<std::size_t>(s)
                                    ].size()
                                ),
                                static_cast<int>(
                                    images[
                                        static_cast<std::size_t>(c)
                                    ].size()
                                ),
                                static_cast<int>(
                                    images[
                                        static_cast<std::size_t>(d)
                                    ].size()
                                )
                            }
                        );
                        for (int cutoff = 0; cutoff <= size; ++cutoff) {
                            const Integer determinant =
                                inner(
                                    images[static_cast<std::size_t>(r)],
                                    images[static_cast<std::size_t>(c)],
                                    cutoff
                                )
                                    * inner(
                                        images[
                                            static_cast<std::size_t>(s)
                                        ],
                                        images[
                                            static_cast<std::size_t>(d)
                                        ],
                                        cutoff
                                    )
                                - inner(
                                      images[
                                          static_cast<std::size_t>(r)
                                      ],
                                      images[
                                          static_cast<std::size_t>(d)
                                      ],
                                      cutoff
                                  )
                                    * inner(
                                        images[
                                            static_cast<std::size_t>(s)
                                        ],
                                        images[
                                            static_cast<std::size_t>(c)
                                        ],
                                        cutoff
                                    );
                            if (
                                !full_initialized
                                || determinant < full_minimum
                            ) {
                                full_initialized = true;
                                full_minimum = determinant;
                                full_r = r;
                                full_s = s;
                                full_c = c;
                                full_d = d;
                                full_cutoff = cutoff;
                            }
                        }
                    }
                }
            }
        }
        std::cout
            << "SU2_LOG_CONCAVE_TAIL_PROFILE"
            << " minimum=" << minimum
            << " q=" << minimum_q
            << " a=" << minimum_a
            << " cutoff=" << minimum_cutoff
            << (minimum < 0 ? " result=COUNTEREXAMPLE"
                            : " result=NO_COUNTEREXAMPLE")
            << " full_minimum=" << full_minimum
            << " full_minor=" << full_r << ',' << full_s
            << ';' << full_c << ',' << full_d
            << " full_cutoff=" << full_cutoff
            << (full_minimum < 0 ? " full_result=COUNTEREXAMPLE"
                                 : " full_result=NO_COUNTEREXAMPLE")
            << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
