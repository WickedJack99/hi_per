#include "dbscan_parallel.h"
#include <chrono>
#include <fstream>
#include <iostream>

using namespace HPC;

int main() {

  std::vector<Point> points = readPointsFromFile("data");

  DBSCAN ds(5, 0.01);
  ds.run(points);

  writePointsToFile(ds.getPoints(), "clustered");

  return 0;
}
