#ifndef DBSCAN_H
#define DBSCAN_H

#include <vector>
#include <set>

#include "point.h"

namespace HPC {

class DBSCAN {
 public:
  DBSCAN(int minPts, double eps);

  void run(const std::vector<Point>& points);

  const std::vector<Point>& getPoints() const { return dataset_; }

 private:
  std::set<int> regionQuery(const Point& point) const;
  bool expandCluster(Point& point, std::set<int>& neighbours, int clusterID);

  // void merge(std::vector<int>& n, const std::vector<int>& nPrime) const;

  const int unclassifiedID = -1;
  const int noiseID = -2;

  const int minPoints_;
  const double epsilon_;

  std::vector<Point> dataset_;
};

}  // namespace HPC

#endif  // DBSCAN_H
