#include <iostream>
#include <fstream>
#include <chrono>
#include "dbscan.h"

using namespace HPC;


int main() {

  std::vector<Point> points = readPointsFromFile("data");

  DBSCAN ds(5, 0.01);
  ds.run(points);

  writePointsToFile(ds.getPoints(), "clustered");

  return 0;
}
