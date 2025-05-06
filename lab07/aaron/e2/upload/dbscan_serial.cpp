#include "dbscan_serial.h"
#include <cmath>
#include <iostream>

namespace HPC {

DBSCAN::DBSCAN(int minPts, double eps) : minPoints_(minPts), epsilon_(eps) {}

void DBSCAN::run(const std::vector<Point>& points) {
  dataset_ = points;
  const int n = dataset_.size();

  int clusterIndex = 0;

  for (int i = 0; i < n; ++i) {
    Point& point = dataset_[i];
    if (point.clusterID < 0) {
      std::set<int> neighbours = regionQuery(point);
      if (neighbours.size() < minPoints_) {
        point.clusterID = noiseID;
      } else {
        clusterIndex++;
        expandCluster(point, neighbours, clusterIndex);
      }
    }
  }

  for (Point point : dataset_) {
    std::cout << point << std::endl;
  }
}

bool DBSCAN::expandCluster(Point& p, std::set<int>& neighbours, int clusterID) {
  p.clusterID = clusterID;

  std::set<int> updatedNeighbours = neighbours;
  neighbours.clear();
  while (updatedNeighbours.size() != neighbours.size()) {
    neighbours = updatedNeighbours;

    for (int i : neighbours) {
      Point& pPrime = dataset_[i];
      if (pPrime.clusterID < 0) {
        pPrime.clusterID = clusterID;  // serves as marking the point as visited
        std::set<int> newNeighbours = regionQuery(pPrime);
        if (newNeighbours.size() >= minPoints_) {
          updatedNeighbours.merge(newNeighbours);
        }
      }
    }
  }
  return true;
}

std::set<int> DBSCAN::regionQuery(const Point& point) const {
  std::set<int> neighbours;
  for (int i = 0; i < dataset_.size(); ++i) {
    if (point.distance(dataset_[i]) <= epsilon_) {
      neighbours.insert(i);
    }
  }
  return neighbours;
}

}  // namespace HPC