// S_h(10^m), m = 1..20, for the parity rock-paper-scissors game G_n(h), in 100-digit decimal arithmetic.
// In G_n(h) a win with the number m pays 2m - h; h = 0 (default) or h = 1 is the first argument.
// S_h(N) is the sum over 3 <= n <= N of the equilibrium probability P(n) of the largest number.
// Complete and partial blocks are summed with the digamma block-sum formula; blocks with
// k >= SERIES_START use a moment expansion with an explicit truncation bound (second column).
// Output: m, S_h(10^m) to 40 decimals, accumulated series error bound.
// Build: g++ -O2 -std=c++17 -pthread alperconstant.cpp -o alperconstant   (Boost headers required)
// Run:   ./alperconstant 0 > S_h0_40dp.txt   and   ./alperconstant 1 > S_h1_40dp.txt
#include <algorithm>
#include <array>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace {

using Real = boost::multiprecision::cpp_dec_float_100;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
const u128 MAX_LIMIT = u128(10'000'000'000'000'000ULL) * 10'000;
constexpr unsigned SERIES_START = 1000;
constexpr unsigned SERIES_ORDER = 8;
unsigned H = 0;  // payment 2m - H

Real to_real(const u128 v) {
    return Real(static_cast<u64>(v >> 64)) * Real("18446744073709551616") + Real(static_cast<u64>(v));
}

void require(const bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error("Check failed: " + description);
}

u128 threshold(const u64 k) {
    // First n of the block k: floor((4k^3+5k)/6)+2 for H = 1, ceil((4k^3+5k)/6)+1 for H = 0.
    const u128 K = k;
    return H == 1 ? (4 * K * K * K + 5 * K) / 6 + 2 : (4 * K * K * K + 5 * K + 5) / 6 + 1;
}

Real block_sum(const u64 k, const u128 n, const std::vector<Real>& central) {
    const u128 first = threshold(k);
    const u128 last = std::min<u128>(n, threshold(k + 1) - 2);
    if (last < first) return 0;

    // This digamma difference sums 1/(2n - H - 2) over the block (the pole with index 1).
    const Real half = Real(H) / 2;
    Real harmonic = (boost::math::digamma(to_real(last) - half)
        - boost::math::digamma(to_real(first) - 1 - half)) / 2;
    Real weighted = 0;
    for (u64 i = 1; i <= k; ++i) {
        const Real residue = 2 * central[i - 1] * (2 * k - 2 * i + 1) * central[k - i];
        weighted += residue * harmonic;
        if (i < k) {
            harmonic += Real(1) / to_real(2 * first - H - 4 * i - 2) + Real(1) / to_real(2 * first - H - 4 * i)
                - Real(1) / to_real(2 * last - H - 4 * i) - Real(1) / to_real(2 * last - H - 4 * i + 2);
        }
    }
    return (to_real(last - first + 1) - Real(2 * k * k + 1) * weighted / 3) / (2 * k + 1);
}

struct Result {
    Real value = 0;
    Real error_bound = 0;
};

std::array<Real, SERIES_ORDER + 1> central_moments(const u64 k) {
    // Moments of 4X/(k-1)-1 for X ~ BetaBinomial(k-1, 1/2, 3/2).
    constexpr int coefficients[SERIES_ORDER + 1][SERIES_ORDER] = {
        {1}, {0}, {1, 2}, {1, 3, 2}, {3, 12, 10, -4},
        {6, 30, 40, 0, -16}, {15, 90, 150, 0, -104, 32},
        {36, 252, 525, 105, -560, -84, 272},
        {91, 728, 1792, 560, -2618, -616, 2248, -544}
    };
    const Real inverse = Real(1) / (k - 1);
    std::array<Real, SERIES_ORDER + 1> moments{};
    moments[0] = 1;
    for (unsigned j = 2; j <= SERIES_ORDER; ++j) {
        for (int degree = static_cast<int>(j) - 1; degree >= 0; --degree) {
            moments[j] = moments[j] * inverse + coefficients[j][degree];
        }
    }
    return moments;
}

Result series_block_sum(const u64 k, const u128 n) {
    const u128 first = threshold(k);
    const u128 last = std::min<u128>(n, threshold(k + 1) - 2);
    if (last < first) return {};
    const Real h = k - 1;
    const Real a = 2 * to_real(first) - H - k - 1;
    const Real b = 2 * to_real(last + 1) - H - k - 1;
    const Real ia = 1 / a, ib = 1 / b;
    const Real ua = h * ia, ub = h * ib;
    const Real ia2 = ia * ia, ib2 = ib * ib;
    const Real ia4 = ia2 * ia2, ib4 = ib2 * ib2;
    const Real ia6 = ia4 * ia2, ib6 = ib4 * ib2;
    const auto moments = central_moments(k);
    Real pa = 1, pb = 1, weighted = 0, bound = 0;
    for (unsigned j = 0; j <= SERIES_ORDER; ++j) {
        const Real integral = j == 0 ? Real(log(b / a) / 2) : Real((pa - pb) / (2 * j));
        const Real second = Real(j + 1) / 6;
        const Real fourth = Real((j + 1) * (j + 2) * (j + 3)) / 90;
        const Real sixth = Real((j + 1) * (j + 2) * (j + 3) * (j + 4) * (j + 5)) / 945;
        const Real remainder = sixth * (pa * ia6 - pb * ib6);
        const Real harmonic = integral + (pa * ia - pb * ib) / 2
            + second * (pa * ia2 - pb * ib2) - fourth * (pa * ia4 - pb * ib4) + remainder;
        weighted += moments[j] * harmonic;
        bound += abs(moments[j]) * remainder;
        pa *= ua;
        pb *= ub;
    }
    const Real rho = 3 * h / a;
    // The last Euler-Maclaurin correction bounds its integral remainder;
    // the reciprocal expansion has a geometric tail bounded by rho^(ORDER+1)/(1-rho).
    bound += to_real(last - first + 1) / a * pow(rho, SERIES_ORDER + 1) / (1 - rho);
    const Real coefficient = Real(2 * k * k + 1) * (2 * k) / (3 * (2 * k + 1));
    return {to_real(last - first + 1) / (2 * k + 1) - coefficient * weighted, coefficient * bound};
}

struct Task {
    u128 n = 0;
    unsigned index = 0;
    unsigned stride = 1;
    const std::vector<Real>* central = nullptr;
    std::vector<Result>* sums = nullptr;
    std::exception_ptr error;
};

void* sum_worker(void* argument) {
    Task& task = *static_cast<Task*>(argument);
    try {
        for (unsigned k = task.index + 1; k < task.sums->size(); k += task.stride) {
            (*task.sums)[k] = k < SERIES_START
                ? Result{block_sum(k, task.n, *task.central), 0}
                : series_block_sum(k, task.n);
        }
    } catch (...) {
        task.error = std::current_exception();
    }
    return nullptr;
}

Result evaluate(const u128 n, unsigned thread_count) {
    require(n <= MAX_LIMIT, "supported limit through 10^20");
    if (n < 3) return {};
    unsigned blocks = 1;
    while (threshold(blocks + 1) <= n) ++blocks;
    std::vector<Real> central(std::min(blocks, SERIES_START));
    std::vector<Result> sums(blocks + 1);
    central[0] = 1;
    for (unsigned j = 1; j < central.size(); ++j) central[j] = central[j - 1] * (2 * j - 1) / (2 * j);
    thread_count = std::min(thread_count, blocks);
    require(thread_count > 0, "positive thread count");
    std::vector<Task> tasks(thread_count);
    std::vector<pthread_t> threads(thread_count);
    unsigned created = 0;
    for (unsigned t = 0; t < thread_count; ++t) {
        tasks[t] = {n, t, thread_count, &central, &sums, {}};
        if (thread_count == 1) {
            sum_worker(&tasks[t]);
        } else {
            if (pthread_create(&threads[t], nullptr, sum_worker, &tasks[t]) != 0) break;
            ++created;
        }
    }
    bool joined = true;
    for (unsigned t = 0; t < created; ++t) joined = pthread_join(threads[t], nullptr) == 0 && joined;
    require(thread_count == 1 || created == thread_count, "pthread_create");
    require(joined, "pthread_join");
    for (const Task& task : tasks) if (task.error) std::rethrow_exception(task.error);
    Result result;
    for (unsigned k = 1; k <= blocks; ++k) {
        result.value += sums[k].value;
        result.error_bound += sums[k].error_bound;
    }
    require(result.error_bound < Real("1e-24"), "series truncation bound");
    return result;
}

unsigned logical_processor_count() {
#ifdef _WIN32
    // All processor groups; hardware_concurrency() counts only the current group (64 of 128 threads here).
    return static_cast<unsigned>(GetActiveProcessorCount(ALL_PROCESSOR_GROUPS));
#else
    return std::thread::hardware_concurrency();
#endif
}

}

int main(int argc, char** argv) {
    H = argc > 1 ? static_cast<unsigned>(std::atoi(argv[1])) : 0;
    require(H <= 1, "h must be 0 or 1");
    u128 n = 1;
    for (unsigned e = 1; e <= 20; ++e) {
        n *= 10;
        const Result r = evaluate(n, logical_processor_count());
        std::cout << e << " " << std::fixed << std::setprecision(40) << r.value << " " << std::scientific << std::setprecision(3) << r.error_bound << "\n" << std::flush;
    }
}
