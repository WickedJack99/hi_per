#include <iostream>
#include <fstream>
#include <chrono>
#include "dbscan.h"

using namespace HPC;

int main()
{

  std::vector<Point> points = readPointsFromFile("data");

  DBSCAN ds(5, 0.01);
  // Zeitmessung starten
  auto start = std::chrono::high_resolution_clock::now();

  ds.run(points);

  // Zeitmessung beenden
  auto end = std::chrono::high_resolution_clock::now();

  // Dauer berechnen in Millisekunden
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout << "Laufzeit: " << duration << " ms" << std::endl;
  writePointsToFile(ds.getPoints(), "clustered");

  return 0;
}
