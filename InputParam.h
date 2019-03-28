#pragma once
#include "basic.h"
#include "threadpool.h"
#include <string>

struct Port
{
	int portID;
	int portType;
	double portCost;
	double portTime;
};

struct legConf
{
    double legCost;
    double legTime;
};

struct Leg
{
	int tailID;
	int headID;
	vector<legConf> confList;
};

struct Demand
{
	int originID;
	int destinationID;
	double weight;
	double revenue;
	double volumeRatio;
	double maxTransTime;
};

struct Scenario
{
	double probability;
	vector<Demand> demandList;
};

class InputParam
{
private:
	string workDir;
public:
	vector<Port> portList;
	vector<Leg> legList;
	vector<Scenario> scenarioList;
	
	double shipWeightCapacity;
	double shipVolumeCapacity;
	int scenarioNum;
	int portNum;
	int legNum;
	int baseStartIndex;
	int baseStartID;
	int baseEndIndex;
	int baseEndID;
	double maxVoyageTime;

	vector<int> fixPortIndices;
	vector<int> addPortIndices;
	vector<pair<int,int>> prioList;
	vector<pair<int,int>> fixPrioList;

	int getPortIndex(int portID);
	int getLegIndex(int tailID, int headID);
	int getPrioIndex(int id_i, int id_j);
	int getFixPrioIndex(int id_i, int id_j);

	string getWorkDir();

	std::threadpool* executor;

	InputParam(string data_dir);

	~InputParam()
	{
		//cout << "d!" << endl;
		delete executor;
	}
};

extern InputParam* param;
