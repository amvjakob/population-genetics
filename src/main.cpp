#include "SimulationsExecutor.hpp"
#include "Data.hpp"
#include <iostream>

int main() {
	
	int n = 400;
	int N = 4000;
	int T = 2000;
	
	/*
	std::cout << "Executing first simulationsExecutor" << std::endl;
	
	SimulationsExecutor simulationsExecutor(n, N, T, {0.4, 0.3, 0.2, 0.1});
	simulationsExecutor.execute();
	* */
	
	std::cout << "Executing second simulationsExecutor" << std::endl;
	
	Data data = Data();
	data.collectAll();
	
	SimulationsExecutor simulationsExecutorFull(data);
	simulationsExecutorFull.execute();
	
	return 0;
}
