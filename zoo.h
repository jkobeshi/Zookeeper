// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#ifndef ZOO_H
#define ZOO_H

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>
#include <queue>
#include <getopt.h>

using namespace std;

class zooData {
public:
	void get_options(int argc, char** argv);
	void read();
	void run();
private:
	enum class Status : uint8_t {
		Wild,
		Safe,
		Border
	};

	struct Cage {
		uint16_t CageID; int32_t x, y; Status area;
		double distance = numeric_limits<double>::infinity(); Cage* edge;
		bool Visited = 0;
	};

	double dist(Cage& a, Cage& b) {
		if ((a.area == Status::Safe && b.area == Status::Wild) ||
			(b.area == Status::Safe && a.area == Status::Wild)) {
			return numeric_limits<double>::infinity();
		}
		else {
			return sqrt(((double(a.x) - double(b.x)) * (double(a.x) - double(b.x))) + ((double(a.y) - double(b.y)) * (double(a.y) - double(b.y))));
		}
	}

	double dist2(Cage& a, Cage& b) {
		return sqrt(((double(a.x) - double(b.x)) * (double(a.x) - double(b.x))) + ((double(a.y) - double(b.y)) * (double(a.y) - double(b.y))));
	}

	double sum() {
		double sum = coords[0].distance;
		for (uint16_t i = 1; i < cord_size; ++i) {
			sum += coords[i].distance;
		} return sum;
	}
	
	bool promising(size_t &permLength) {
		if ((Path.size() - permLength) <= 4) {
			return true;
		}
		if (cost >= bestCost) {
			return false;
		}
		for (uint16_t i = uint16_t(permLength + 1); i < cord_size; ++i) {
			coords[Path[i]].distance = numeric_limits<double>::infinity();
			coords[Path[i]].Visited = 0;
		}
		coords[Path[0]].distance = dMatrix[Path[permLength]][Path[0]];
		coords[Path[permLength - 1]].distance = dMatrix[Path[permLength]][Path[permLength - 1]];
		for (size_t i = permLength; i < Path.size(); ++i) {
			if (dMatrix[Path[i]][Path[0]] < coords[Path[0]].distance) {
				coords[Path[0]].distance = dMatrix[Path[i]][Path[0]];
			}
			if (dMatrix[Path[i]][Path[permLength - 1]] < coords[Path[permLength - 1]].distance) {
				coords[Path[permLength - 1]].distance = dMatrix[Path[i]][Path[permLength - 1]];
			}
		}
		size_t prev = permLength, min = permLength; double MST = 0.0;
		for (size_t nVisit = permLength + 1; nVisit < Path.size(); ++nVisit) {
			for (size_t i = permLength + 1; i < Path.size(); ++i) {
				if (!coords[Path[i]].Visited) {
					if (dMatrix[Path[i]][Path[prev]] < coords[Path[i]].distance) {
						coords[Path[i]].distance = dMatrix[Path[i]][Path[prev]];
					}
					if (coords[Path[min]].Visited || (coords[Path[min]].distance > coords[Path[i]].distance)) {
						min = i;
					}
				}
			}
			prev = min, coords[Path[prev]].Visited = 1; MST += coords[Path[prev]].distance;
		}
		return ((cost + coords[Path[0]].distance + coords[Path[permLength - 1]].distance + MST) < bestCost);
	}


	void genPerms(size_t permLength) {
		if (permLength == Path.size()) {
			cost += dMatrix[Path[permLength]][Path[permLength - 1]];
			if (cost < bestCost) {
				bestCost = cost;
				bestPath = Path;
			}
			cost -= dMatrix[Path[permLength]][Path[permLength - 1]];
			return;
		} // if
		if (!promising(permLength))
			return;
		for (size_t i = permLength; i < Path.size(); ++i) {
			swap(Path[permLength], Path[i]);
			cost += dMatrix[Path[permLength]][Path[permLength - 1]];
			genPerms(permLength + 1);
			cost -= dMatrix[Path[permLength]][Path[permLength - 1]];
			swap(Path[permLength], Path[i]);
		} // for
	} // genPerms()

	void MST();
	void MSTprint();
	void FASTTSP();
	void OPTTSP();
	void TSPprint();
	string format;
	uint16_t safe = 0, border = 0, wild = 0, cord_size = 0;
	vector<Cage> coords;
	vector<vector<double>> dMatrix;
	//OPTTSP
	vector<uint16_t> Path; vector<uint16_t> bestPath;
	double cost = 0.0; double bestCost = 0.0;
};

void zooData::get_options(int argc, char** argv) {
	int gotopt;
	int option_index = 0;
	option long_opts[] = {
		{ "mode", required_argument, nullptr, 'm'},
		{ "help", no_argument, nullptr, 'h'} };
	while ((gotopt = getopt_long(argc, argv, "m:h", long_opts, &option_index)) != -1) {
		switch (gotopt) {
			case 'm':
				format = optarg;
				break;
			case 'h':
				cout << "Use the command -m or --mode to execute the desired mode!\n";
				exit(0);
				break;
			default:
				cerr << "Invalid command line option\n";
				exit(1);
				break;
		}
	}
	if (!(format == "MST" || format == "OPTTSP" || format == "FASTTSP")) {
		cerr << "Invalid mode\n";
		exit(1);
	}
}

void zooData::read() {
	cin >> cord_size;
	coords.resize(cord_size);
	uint16_t i = 0;
	if (format == "MST") {
		while (cin >> coords[i].x >> coords[i].y) {
			coords[i].CageID = i;
			if (coords[i].x <= 0 && coords[i].y <= 0) {
				if (coords[i].x == 0 || coords[i].y == 0) {
					coords[i].area = Status::Border;
					++border;
				}
				else {
					coords[i].area = Status::Wild;
					++wild;
				}
			} 
			else {
				coords[i].area = Status::Safe;
			}++i;
		}
	}
	else {
		while (cin >> coords[i].x >> coords[i].y) {
			coords[i].CageID = i;
			++i;
		}
	}
}

void zooData::run() {
	if (format == "MST") {
		if (safe > 0 && wild > 0 && border == 0) {
			cerr << "Cannot construct MST\n";
			exit(1);
		}
		MST();
		MSTprint();
	}
	else if (format == "OPTTSP") {
		OPTTSP();
	}
	else {
		FASTTSP();
		TSPprint();
	}
}

#endif