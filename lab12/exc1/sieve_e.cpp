#include <vector>
#include <iostream>
#include <mpi.h>
#include <limits>
#include <fstream>

std::vector<int> get_start_condition(int N, int count_proc, int rank)
{
  int total_elements = N - 1;
  int base_chunk = total_elements / count_proc;
  int remainder = total_elements % count_proc;

  int start_offset = rank * base_chunk + std::min(rank, remainder);
  int chunk_size = base_chunk + (rank < remainder ? 1 : 0);

  int start = 2 + start_offset;
  int end = start + chunk_size;

  std::vector<int> start_set;
  for (int i = start; i < end; ++i)
  {
    if (i == 2 || i % 2 == 1)
      start_set.emplace_back(i);
  }

  return start_set;
}

// mark element as not prim by setting to 0
// point out next k while sieving
int sieve(std::vector<int> &set, int k)
{
  bool next_k_set = false;
  int next_k = k;
  for (int &elem : set)
  {
    if ((elem > k) && (elem != 0))
    {
      if (!next_k_set)
      {
        next_k = elem;
        next_k_set = true;
      }
      if (elem % k == 0)
      {
        elem = 0;
      }
    }
  }
  return next_k;
}

void print_vec(std::vector<int> vec)
{
  for (int elem : vec)
  {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;
}

void save_result_to_file(std::vector<int> vec, int N)
{
  std::ofstream outfile("primes_" + std::to_string(N) + ".txt");
  for (int val : vec)
  {
    outfile << val << ",";
  }
  outfile << std::endl;
  outfile.close();
}

void save_runtime_to_file(double runtime, int N, int numProc)
{
  std::ofstream outfile("runtime_p_" + std::to_string(numProc) + "_N_" + std::to_string(N) + ".txt");
  outfile << runtime << "," << N;
  outfile << std::endl;
  outfile.close();
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  double start_time = MPI_Wtime();
  int N = atoi(argv[1]);
  int numProc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> list_to_sieve = get_start_condition(N, numProc, rank);

  int k = 2;
  int global_min_k = 0;
  int local_max_N = list_to_sieve[list_to_sieve.size() - 1];
  while (true)
  {
    int offered_k = (k < local_max_N) ? sieve(list_to_sieve, k) : std::numeric_limits<int>::max();

    MPI_Allreduce(&offered_k, &global_min_k, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (global_min_k * global_min_k >= N)
      break;

    k = global_min_k;
  }

  std::vector<int> recvcounts(numProc);
  std::vector<int> displs(numProc);

  int local_size = list_to_sieve.size();
  std::vector<int> final_results(numProc * local_size);

  // 1. Share local sizes from all ranks to rank 0
  MPI_Gather(&local_size, 1, MPI_INT,
             recvcounts.data(), 1, MPI_INT,
             0, MPI_COMM_WORLD);

  // 2. Calculate displacements on rank 0
  if (rank == 0)
  {
    displs[0] = 0;
    for (int i = 1; i < numProc; ++i)
    {
      displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    // Resize final_results based on total size
    int total = displs[numProc - 1] + recvcounts[numProc - 1];
    final_results.resize(total);
  }

  if (rank == 0)
  {
    std::copy(list_to_sieve.begin(), list_to_sieve.end(), final_results.begin());
    MPI_Gatherv(MPI_IN_PLACE, local_size, MPI_INT,
                final_results.data(), recvcounts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);
  }

  else
  {
    MPI_Gatherv(list_to_sieve.data(), local_size, MPI_INT,
                nullptr, nullptr, nullptr, MPI_INT,
                0, MPI_COMM_WORLD);
  }

  double end_time = MPI_Wtime();
  double elapsed = end_time - start_time;

  if (rank == 0)
  {
    std::vector<int> non_zero_results;
    for (int val : final_results)
    {
      if (val != 0)
        non_zero_results.push_back(val);
    }

    save_result_to_file(non_zero_results, N);
    save_runtime_to_file(elapsed, N, numProc);
  }

  MPI_Finalize();
  return 0;
}