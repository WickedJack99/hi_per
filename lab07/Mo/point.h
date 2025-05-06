#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>

/**
 * Class representing a point in the dataset.
 *
 * Stores the coordinates of the point, its cluster ID, and whether it is a core
 * point.
 */
class Point {
 public:
  Point(const std::vector<double>& coordinatesIn);

  double& operator()(int i);
  const double& operator()(int i) const;

  double distance(const Point& other) const;

  std::vector<double> coordinates;
  int clusterID = -1;
  bool isCorePoint = false;
};

/**
 * Read points from a file and return them as a vector of Point objects.
 */
std::vector<Point> readPointsFromFile(const std::string& filename);

/**
 * Print a point to an output stream. The 
 * coordinates are separated by tabs, and the
 * cluster ID is printed at the end.
 */
std::ostream& operator<<(std::ostream& os, const Point& point);

/**
 * Write points to a file.
 * 
 * Each point is written on a new line, with
 * coordinates separated by tabs and the
 * cluster ID at the end.
 * 
 * Can be read with numpy.loadtxt, the last column give the cluster ID.
 */
void writePointsToFile(const std::vector<Point>& points,
                       const std::string& filename);

#endif  // POINT_H