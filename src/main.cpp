#include "SimulationsExecutor.hpp"
#include "Data.hpp"

int main() {
	/*
	int n = 10;
	int N = 1000;
	int T = 1000;
	
	SimulationsExecutor simulationsExecutor(n, N, T, {0.4, 0.3, 0.2, 0.1});
	simulationsExecutor.execute();
	*/
	
	Data data = Data();
	data.collectAll();
	
	SimulationsExecutor simulationsExecutorFull(data);
	simulationsExecutorFull.execute();
	
	
	return 0;
}
