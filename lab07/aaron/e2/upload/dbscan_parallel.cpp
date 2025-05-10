#include "dbscan_parallel.h"
#include <atomic>
#include <cmath>
#include <iostream>
#include <omp.h>

namespace HPC {

DBSCAN::DBSCAN(int minPts, double eps) : minPoints_(minPts), epsilon_(eps) {}

void DBSCAN::run(const std::vector<Point> &points, int threadNum) {

  dataset_ = points;
  const int n = dataset_.size();

  initializeNeighbors(threadNum);

  int clusterIndex = 0;
  for (int i = 0; i < n; ++i) {
    Point &point = dataset_[i];
    if (point.clusterID < 0) {
      std::set<int> neighbours = point.neighbors;
      if (neighbours.size() < minPoints_) {
        point.clusterID = noiseID;
      } else {
        clusterIndex++;
        expandCluster(point, neighbours, clusterIndex);
      }
    }
  }
}

bool DBSCAN::expandCluster(Point &p, std::set<int> &neighbours, int clusterID) {
  p.clusterID = clusterID;

  std::set<int> updatedNeighbours = neighbours;

  // Use of do-while instead of clearing neighbors
  do {
    neighbours = updatedNeighbours;

    for (int i : neighbours) {
      Point &pPrime = dataset_[i];
      if (pPrime.clusterID < 0) {
        pPrime.clusterID = clusterID; // serves as marking the point as visited
        std::set<int> newNeighbours = pPrime.neighbors;
        if (newNeighbours.size() >= minPoints_) {
          updatedNeighbours.merge(newNeighbours);
        }
      }
    }
  } while (updatedNeighbours.size() != neighbours.size());
  return true;
}

void DBSCAN::initializeNeighbors(int threadNum) {
#pragma omp parallel for num_threads(threadNum)
  for (int i = 0; i < dataset_.size(); ++i) {
    Point &pointToCheckNeighborsFor = dataset_[i];
    for (int j = 0; j < dataset_.size(); ++j) {
      if (pointToCheckNeighborsFor.distance(dataset_[j]) <= epsilon_) {
        pointToCheckNeighborsFor.neighbors.insert(j);
      }
    }
  }
}

} // namespace HPC
