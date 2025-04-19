#include <iostream>
#include <thread>
#include <mutex>

std::mutex fooPrinted;
std::mutex barPrinted;

const int n = 100;
void foo() {
  for (int i = 0; i < n; ++i) {
    barPrinted.lock();
    std::cout << "foo";
    fooPrinted.unlock();
  }
}
void bar() {
  for (int i = 0; i < n; ++i) {
    fooPrinted.lock();
    std::cout << "bar\n";
    barPrinted.unlock();
  }
}
int main() {
  fooPrinted.lock();
  std::thread t1(foo);
  std::thread t2(bar);
  t1.join();
  t2.join();
}