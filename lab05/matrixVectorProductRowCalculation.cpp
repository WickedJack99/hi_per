#include <benchmark/benchmark.h>
#include <cmath>
#include <iostream>
#include <thread>
#include "matrix.h"
#include "test.h"
#include <fstream>
#include <string>

// Create the matrix and vector to be multiplied and fill them
// with some sensible initial values.
std::pair<Matrix, std::vector<double>> createMatrixAndVector()
{
    const int n = 1e3 * 9;
    Matrix mat(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            mat(i, j) = pow(-1, i) * (i + j);
        }
    }

    std::vector<double> vec(n);
    for (int i = 0; i < n; ++i)
    {
        vec[i] = 1. / (i + 1);
    }

    return std::pair(mat, vec);
}

// Verify that the computed result is correct. Rather inefficient,
// since it runs on a single core.
void verifyResult(const std::vector<double> result)
{
    auto [mat, vec] = createMatrixAndVector();
    const int n = vec.size();

    for (int i = 0; i < n; ++i)
    {
        double expected = 0;
        for (int j = 0; j < n; ++j)
        {
            expected += mat(i, j) * vec[j];
        }
        check(result[i], expected);
    }
}

void computeResult(
    const Matrix &mat,
    const std::vector<double> &vec,
    std::vector<double> &result,
    int row)
{
    int n = vec.size();
    for (int col = 0; col < n; col++)
    {
        result[row] += mat(row, col) * vec[col];
    }
}

void workerThread(
    const Matrix &mat,
    const std::vector<double> &vec,
    std::vector<double> &result,
    std::atomic<int> &rowIndex,
    int &rowsComputed)
{
    while (true)
    {
        int row = rowIndex.fetch_add(1);
        if (row >= vec.size()) break;
        computeResult(mat, vec, result, row);
        rowsComputed++;
    }
}

void storeRowsComputedPerThread(std::vector<int> &rowsComputedPerThread, int threadCount, int stateCount) {
    std::stringstream fileName;
    fileName << "stats" << threadCount << "_" << stateCount << ".csv";
    std::ofstream MyFile(fileName.str());
    for (auto entry : rowsComputedPerThread) {
      MyFile << entry << ",";
    }
    MyFile.close();
}

void benchmarkComputeResult(benchmark::State &state)
{
    int threadCount = state.range(0);
    auto [mat, vec] = createMatrixAndVector();
    int stateCount = 0;
    for (auto _ : state)
    {
        std::atomic<int> rowIndex = 0;
        std::vector<std::thread> threadPool;
        std::vector<double> result(vec.size(), 0);
        std::vector<int> rowsComputedPerThread(threadCount, 0);
        for (int threadIndex = 0; threadIndex < threadCount; threadIndex++)
        {
            threadPool.emplace_back(
                workerThread, std::ref(mat), std::ref(vec),
                std::ref(result), std::ref(rowIndex), std::ref(rowsComputedPerThread[threadIndex]));
        }

        for (int threadIndex = 0; threadIndex < threadCount; threadIndex++)
        {
            threadPool[threadIndex].join();
        }
        storeRowsComputedPerThread(rowsComputedPerThread, threadCount, stateCount);
        benchmark::DoNotOptimize(result);
        stateCount++;
    }
}

int main(int argc, char **argv)
{
    ::benchmark::Initialize(&argc, argv);

    std::vector<int> threads = {1, 2, 3, 4, 6, 9, 12};

    for (int i = 0; i < threads.size(); i++)
    {
        benchmark::RegisterBenchmark("idk", benchmarkComputeResult)
            ->Arg(threads[i])
            ->Unit(benchmark::kMillisecond);
    }

    ::benchmark::RunSpecifiedBenchmarks();

    return 0;
}