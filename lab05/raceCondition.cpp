// adapted from modernescpp.com

#include <iostream>
#include <thread>
#include "test.h"
#include <atomic>
#include <mutex>

std::mutex mut;

struct Account {
  std::atomic<int> balance{100};
};

void transferMoney(int amount, Account& from, Account& to) {
  using namespace std::chrono_literals;
  mut.lock();
  if (from.balance.load() >= amount) {
    std::this_thread::sleep_for(1ns);
    from.balance -= amount;
    std::this_thread::sleep_for(1ns);
    to.balance += amount;
  }
  mut.unlock();
}

int testTransferMoney() {
  Account account1;
  Account account2;

  std::thread thr1(transferMoney, 80, std::ref(account1), std::ref(account2));
  std::thread thr2(transferMoney, 60, std::ref(account1), std::ref(account2));
  std::thread thr3(transferMoney, 10, std::ref(account2), std::ref(account1));

  thr1.join();
  thr2.join();
  thr3.join();

  std::cout << "\nChecking balance of account1: ";
  check(account1.balance > 0, true);
  std::cout << "Checking balance of account2: ";
  check(account2.balance > 0, true);
  std::cout << "Checking sum of accounts:     ";
  check(account1.balance + account2.balance, 200);

  return account1.balance + account2.balance;
}

int main() {
  int erroneousTransfers = 0;
  for (int i = 0; i < 10000; ++i) {
    int sum = testTransferMoney();
    if (sum != 200) {
      erroneousTransfers++;
    }
  }

  std::cout << "\n\nThere were " << erroneousTransfers
            << " transfers where money appeared or disappeared.\n";
}