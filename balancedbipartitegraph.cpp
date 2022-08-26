#include "balancedbipartitegraph.h"

BalancedBipartiteGraph::BalancedBipartiteGraph(int n_): n{n_} {}

BalancedBipartiteGraph::BalancedBipartiteGraph(int n_, std::vector<std::pair<int,int>> edges_) {
  this->n = n_;
  //this->edges = edges_; //VORSICHT!!!! HIER WIRD NUR DIE ADRESSE KOPIERT, NICHT DER INHALT???????
  for(auto e : edges_) {
	this->edges.push_back(e);
  }
}

bool BalancedBipartiteGraph::findMatch(int u, std::vector<bool>& visited, std::vector<int>& matchL, std::vector<int>& matchR) {
  for(int v = 0; v <= this->n - 1; v++) {
    std::vector<std::pair<int,int>>::iterator it = std::find(this->edges.begin(), this->edges.end(), std::make_pair(u,v));

    if((it != this->edges.end()) && !visited[v]) {
      visited[v] = true;

      if(matchR[v] < 0 || findMatch(matchR[v], visited, matchL, matchR)) {
        matchL[u] = v;
        matchR[v] = u;

		//std::cout << "matchL: " << matchL[0] << ", " << matchL[1] << ", " << matchL[2] << std::endl;
		//std::cout << "matchR: " << matchR[0] << ", " << matchR[1] << ", " << matchR[2] << std::endl;

		return true;
      }
    }
  }

  //std::cout << "visited: " << visited[0] << ", " << visited[1] << ", " << visited[2] << std::endl;
  
  return false;
}


std::vector<int> BalancedBipartiteGraph::findMaxMatch() {
	std::vector<bool> visited(this->n, false);
	std::vector<int> mL(this->n, -1);
	std::vector<int> mR(this->n, -1);


	for(int u = 0; u <= this->n - 1; u++) {
	for(int v = 0; v <= this->n - 1; v++) {
	  visited[v] = false;
	}

	this->findMatch(u, visited, mL, mR);
	}

	for(auto x : mL) {
	std::cout << "mL: " << x << std::endl;
	}

	return mL;
}