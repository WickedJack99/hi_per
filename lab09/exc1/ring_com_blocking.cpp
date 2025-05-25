#include <iostream>
#include <mpi.h>
#include <thread>
#include <chrono>

void printSentValue(int rank, int receiverId) {
  std::cout << "Proc " << rank << " sent to " << receiverId << std::endl;
}

void printReceiveValue(int rank, int senderId) {
  std::cout << "Proc " << rank << " received from " << senderId << std::endl;
}

int ringSendRecv(int sendValue) {
  int numProc;
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  int receiverId = (sendValue + 1) % numProc;
  int senderId = (sendValue - 1) % numProc;
  int recvValue;
  if (sendValue == 0) {
    MPI_Send(&sendValue, 1, MPI_INT, receiverId, 0, MPI_COMM_WORLD);
    printSentValue(sendValue, receiverId);
    MPI_Recv(&recvValue, 1, MPI_INT, senderId, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Sleep not necessary, only for pretty .out file
    std::this_thread::sleep_for(std::chrono::seconds(1)); 
    printReceiveValue(sendValue, recvValue);
  } else {
    MPI_Recv(&recvValue, 1, MPI_INT, senderId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::this_thread::sleep_for(std::chrono::seconds(1)); 
    printReceiveValue(sendValue, recvValue);
    MPI_Send(&sendValue, 1, MPI_INT, receiverId, 0, MPI_COMM_WORLD);
    printSentValue(sendValue, receiverId);
  }
  return recvValue;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ringSendRecv(rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  std::cout << "Process on " << processor_name << ", " << rank << "/" << numProc << " finished.\n";
  MPI_Finalize();
}