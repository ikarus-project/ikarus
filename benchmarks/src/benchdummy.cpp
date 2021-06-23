#include <benchmark/benchmark.h>
#include <ranges>
#include <vector>
bool isPrime(int n) {
  if (n <= 1) return false;

  for (int i = 2; i < n; i++)
    if (n % i == 0) return false;

  return true;
}

static void ranges1(benchmark::State& state) {
  auto v = std::ranges::iota_view{1, state.range(0)};

  for (auto _ : state) {
    for (int i : v | std::views::filter(isPrime)) benchmark::DoNotOptimize(i);
  }
}
// Register the function as a benchmark
BENCHMARK(ranges1)->RangeMultiplier(2)->Range(8, 8 << 10);

static void ranges2(benchmark::State& state) {
  auto v = std::ranges::iota_view{1, state.range(0)};

  for (auto _ : state) {
    for (int i : std::ranges::filter_view(v, isPrime)) benchmark::DoNotOptimize(i);
  }
}
// Register the function as a benchmark
BENCHMARK(ranges2)->RangeMultiplier(2)->Range(8, 8 << 10);

static void byhand(benchmark::State& state) {
  auto v = std::ranges::iota_view{1, state.range(0)};

  for (auto _ : state) {
    for (int i : v)
      if (isPrime(i)) benchmark::DoNotOptimize(i);
  }
}
BENCHMARK(byhand)->RangeMultiplier(2)->Range(8, 8 << 10);
