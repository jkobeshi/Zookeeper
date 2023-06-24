// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#include "zoo.h"

void zooData::MST() {
	coords[0].Visited = 1, coords[0].distance = 0.0; uint16_t prev = 0, min = 0;
	for(uint16_t nVisit = 1; nVisit < cord_size; ++nVisit) {
		for (uint16_t i = 1; i < cord_size; ++i) {
			if (!coords[i].Visited) { 
				if (dist(coords[i], coords[prev]) < coords[i].distance) {
					coords[i].distance = dist(coords[i], coords[prev]);
					coords[i].edge = &coords[prev];
				}
				if (coords[min].Visited || (coords[min].distance > coords[i].distance)) {
					min = i;
				}
			}
		} prev = min, coords[prev].Visited = 1;
	}
}

void zooData::MSTprint() {
	cout << setprecision(2) << fixed << sum() << "\n";
	for (uint16_t i = 1; i < cord_size; ++i) {
		if (i < coords[i].edge->CageID) {
			cout << i << " " << coords[i].edge->CageID << "\n";
		}
		else {
			cout << coords[i].edge->CageID << " " << i << "\n";
		}
	}
}

void zooData::FASTTSP() {
	coords[0].Visited = 1; coords[0].edge = &coords[0];
	uint16_t prev = 0, insert;
	for (uint16_t i = 1; i < cord_size; ++i) {
		prev = 0; insert = coords[0].edge->CageID;
		while (insert != 0) {
			if (((dist2(coords[insert], coords[i]) + dist2(*coords[insert].edge, coords[i])) - coords[insert].distance)
				< ((dist2(coords[prev], coords[i]) + dist2(*coords[prev].edge, coords[i])) - coords[prev].distance)) {
				prev = insert;
			}
			insert = coords[insert].edge->CageID;
		}
		coords[i].edge = coords[prev].edge, coords[prev].edge = &coords[i];
		coords[i].distance = dist2(coords[i], *coords[i].edge);
		coords[prev].distance = dist2(coords[prev], *coords[prev].edge);
		coords[i].Visited; prev = i;
	}
}

void zooData::OPTTSP() {
	dMatrix.resize(cord_size, vector<double>(cord_size));
	for (uint16_t i = 0; i < cord_size; ++i) {
		for (uint16_t j = 0; j < cord_size; ++j) {
			dMatrix[i][j] = dist2(coords[i], coords[j]);
		}
		Path.push_back(i);
	}
	coords[0].Visited = 1; coords[0].edge = &coords[0];
	uint16_t prev = 0, insert;
	for (uint16_t i = 1; i < cord_size; ++i) {
		prev = 0; insert = coords[0].edge->CageID;
		while (insert != 0) {
			if (((dMatrix[insert][i] + dMatrix[coords[insert].edge->CageID][i]) - coords[insert].distance)
				< ((dMatrix[prev][i] + dMatrix[coords[prev].edge->CageID][i]) - coords[prev].distance)) {
				prev = insert;
			}
			insert = coords[insert].edge->CageID;
		}
		coords[i].edge = coords[prev].edge, coords[prev].edge = &coords[i];
		coords[i].distance = dMatrix[i][coords[i].edge->CageID];
		coords[prev].distance = dMatrix[prev][coords[prev].edge->CageID];
		coords[i].Visited; prev = i;
	}
	
	int16_t path = coords[0].edge->CageID; 
	bestPath.push_back(0), bestCost = sum();
	while (path != 0) {
		bestPath.push_back(path);
		path = coords[path].edge->CageID;
	} 
	genPerms(1);
	cout << setprecision(2) << fixed << bestCost << "\n0";
	for (uint16_t i = 1; i < cord_size; ++i) {
		cout << " " << bestPath[i];
	} cout << endl;
}

void zooData::TSPprint() {
	cout << setprecision(2) << fixed << sum() << "\n0";
	uint16_t path = 0;
	while (coords[path].edge->CageID != 0) {
		cout << " " << coords[path].edge->CageID;
		path = coords[path].edge->CageID;
	}
	cout << "\n";
}

int main(int argc, char* argv[]) {
	zooData zoo;
	zoo.get_options(argc, argv);
	zoo.read();
	zoo.run();
}