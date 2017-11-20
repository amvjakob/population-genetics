#include "SimulationsExecutor.hpp"
#include "Random.hpp"

#include <thread>
#include <iostream>

void hello(int i) {
	
	std::cout << "Hello world my guy " << i << std::endl;
}


int main() {
	
	/*
	std::thread threads[100];
	
	for (int i = 0; i < 100; i++) {
		threads[i] = std::thread(hello, i);
	}
	
	for (int i = 0; i < 100; i++) {
		threads[i].join();
	}
	* */

	
	int n = 500;
	int N = 5000;
	int T = 3000;
	
	
	SimulationsExecutor simulationsExecutor(n, N, T, {0.4, 0.3, 0.3});
	simulationsExecutor.execute();
	
	
	return 0;
}
