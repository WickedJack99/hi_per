#include <mpi.h>
#include <omp.h>
#include <vector>
#include <fstream>

int main(int argc, char *argv[])
{
    int required = MPI_THREAD_MULTIPLE;
    int provided;
    MPI_Init_thread(&argc, &argv, required, &provided);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    int leftNeighbor = (rank - 1) % numProc;
    int rightNeighbor = (rank + 1) % numProc;

    std::ofstream outfile("rank_" + std::to_string(rank) + ".txt");

#pragma omp parallel
    {
        int threadId = omp_get_thread_num();

        std::vector<int> rankThreadId = {rank, threadId};
        std::vector<int> recvBuffer(2);

        MPI_Request sendReq, recvReq;
        MPI_Status statusRec, statusSend;
        MPI_Isend(rankThreadId.data(), 2, MPI_INT, rightNeighbor, threadId, MPI_COMM_WORLD, &sendReq);
        MPI_Irecv(recvBuffer.data(), 2, MPI_INT, leftNeighbor, threadId, MPI_COMM_WORLD, &recvReq);
        MPI_Wait(&sendReq, &statusSend);

#pragma omp critical
        outfile << " sent (" << rankThreadId[0] << "," << rankThreadId[1] << ")" << std::endl;

        MPI_Wait(&recvReq, &statusRec);
#pragma omp critical
        outfile << " received (" << recvBuffer[0] << "," << recvBuffer[1] << ")" << std::endl;
    }

    outfile.close();

    MPI_Finalize();
}