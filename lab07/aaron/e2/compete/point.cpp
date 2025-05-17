#include <fstream>
#include <iostream>

#include "point.h"

Point::Point(const std::vector<double> &coordinatesIn)
    : coordinates(coordinatesIn) {}

double &Point::operator()(int i) { return coordinates[i]; }

const double &Point::operator()(int i) const { return coordinates[i]; }

double Point::distance(const Point &other) const {
  double distance = 0;
  for (int i = 0; i < coordinates.size(); ++i) {
    const double p = coordinates[i];
    const double q = other.coordinates[i];
    distance += (p - q) * (p - q);
  }

  return distance;
}

std::vector<Point> readPointsFromFile(const std::string &filename) {
  std::vector<Point> points;
  std::ifstream fin(filename);

  double x, y;

  while (fin >> x >> y) {
    Point point({x, y});
    points.push_back(point);
  }
  return points;
}

std::ostream &operator<<(std::ostream &os, const Point &point) {
  for (auto coordinate : point.coordinates) {
    os << coordinate << "\t";
  }
  os << point.clusterID;
  os << "\t" << point.neighbors.size();
  return os;
}

void writePointsToFile(const std::vector<Point> &points,
                       const std::string &filename) {
  std::ofstream fout(filename);
  for (auto point : points) {
    fout << point << "\n";
  }
}
