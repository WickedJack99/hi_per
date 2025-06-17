#include <vector>
#include <iostream>
#include <mpi.h>

std::vector<int> get_start_condition(int N, int count_proc, int rank)
{
    int total_elements = N - 1; // because range starts at 2
    int base_chunk = total_elements / count_proc;
    int remainder = total_elements % count_proc;

    int start_offset = rank * base_chunk + std::min(rank, remainder);
    int chunk_size = base_chunk + (rank < remainder ? 1 : 0);

    int start = 2 + start_offset;
    int end = start + chunk_size; // exclusive

    std::vector<int> start_set;
    for (int i = start; i < end; ++i)
    {
        start_set.emplace_back(i);
    }

    return start_set;
}


// mark element as not prim by setting to 0

// 2,3,4
// 0,1,2
// even numbers are at even indices
// index is (number - 2) or (number - start)

void sieve(std::vector<int> &set, int k)
{
  for (int &elem : set)
  {
    if ((elem != k) && (elem != 0) && (elem % k == 0))
    {
      elem = 0;
    }
  }
}

void print_vec(std::vector<int> vec)
{
  for (int elem : vec)
  {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;
}

int main(int argc, char *argv[])
{
  // MPI_Init(&argc, &argv);
  int procs = atoi(argv[1]);
  int N = atoi(argv[2]);
  std::cout << "Procs: " << procs << " N: " << N << std::endl;
  for (int proc = 0; proc < procs; proc++)
  {
    std::vector<int> elems = get_start_condition(N, procs, proc);
    print_vec(elems);
  }

  // MPI_Finalize();
}