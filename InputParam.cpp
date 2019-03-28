#include "InputParam.h"
#include "ReadCsv.h"

InputParam *param = NULL;

int InputParam::getPortIndex(int portID)
{
    int index = -1;
    for (int i = 0; i < portList.size(); i++)
    {
        if (portID == portList[i].portID)
        {
            index = i;
            break;
        }
    }
    return index;
}

int InputParam::getLegIndex(int tailID, int headID)
{
    int index = -1;
    for (int i = 0; i < legList.size(); i++)
    {
        if ((tailID == legList[i].tailID) && (headID == legList[i].headID))
        {
            index = i;
            break;
        }
    }
    return index;
}

int InputParam::getPrioIndex(int id_i, int id_j)
{
    int index = -1;
    for (int i = 0; i < prioList.size(); i++)
    {
        if (prioList[i].first == id_i && prioList[i].second == id_j)
        {
            index = i;
            break;
        }
    }
    return index;
}

int InputParam::getFixPrioIndex(int id_i, int id_j)
{
    int index = -1;
    for (int i = 0; i < fixPrioList.size(); i++)
    {
        if (fixPrioList[i].first == id_i && fixPrioList[i].second == id_j)
        {
            index = i;
            break;
        }
    }
    return index;
}

string InputParam::getWorkDir()
{
    return workDir;
}

InputParam::InputParam(string data_dir)
{
    workDir = data_dir;
    //read data
    vector<vector<double>> portData = ReadCsv::read(workDir + "port.csv");
    vector<vector<double>> legData = ReadCsv::read(workDir + "leg.csv");
    vector<vector<double>> demandData = ReadCsv::read(workDir + "demand.csv");
    vector<vector<double>> probData = ReadCsv::read(workDir + "probability.csv");
    vector<vector<double>> paramData = ReadCsv::read(workDir + "param.csv");

	if (portData.size() == 0 || legData.size() == 0 || demandData.size() == 0 || paramData.size() == 0)
		throw("no data record!");

    //initialize params
	shipWeightCapacity = paramData[0][0];
	shipVolumeCapacity = paramData[1][0];
	maxVoyageTime = paramData[2][0];

    portNum = portData.size();
    legNum = legData.size();

    portList = vector<Port>(portNum);
    legList = vector<Leg>(legNum);
    scenarioList = vector<Scenario>();

    fixPortIndices = vector<int>();
    addPortIndices = vector<int>();
    prioList = vector<pair<int, int>>();
    fixPrioList = vector<pair<int, int>>();

    // port list
    for (int i = 0; i < portNum; i++)
    {
        Port p;
        p.portID = Round(portData[i][0]);
        p.portType = Round(portData[i][1]);
        p.portCost = portData[i][2];
        p.portTime = portData[i][3];
        portList[i] = p;
        if (p.portType == 0)
        {
            baseStartIndex = i;
            baseStartID = p.portID;
        }
        if (p.portType == -1)
        {
            baseEndIndex = i;
            baseEndID = p.portID;
        }
        if (p.portType == 1)
            fixPortIndices.push_back(i);
        if (p.portType == 2)
            addPortIndices.push_back(i);
    }

    for (int i = 0; i < portNum; i++)
    {
        if (portList[i].portType <= 0) // base port
            continue;
        for (int j = 0; j < portNum; j++)
        {
            if (i == j)
                continue;
            if (portList[j].portType <= 0)
                continue;
            if ((portList[i].portType == 1) && (portList[j].portType == 1))
            {
                fixPrioList.emplace_back(portList[i].portID, portList[j].portID);
            }
            prioList.emplace_back(portList[i].portID, portList[j].portID);
        }
    }

    // leg list
    for (int i = 0; i < legData.size(); i++)
    {
        Leg l;
        l.tailID = Round(legData[i][0]);
        l.headID = Round(legData[i][1]);
        l.confList = vector<legConf>();
        int confNum = (legData[i].size() - 2)/2;
        for (int j = 0; j < confNum; j++)
        {
            double legCost = legData[i][2*j + 2];
            double legTime = legData[i][2*j + 3];
            if (legCost < 0)
                continue;
            legConf lc;
            lc.legCost = legCost;
            lc.legTime = legTime;
            l.confList.push_back(lc);
        }
        legList[i] = l;
    }

    // scenario list
    for (int i = 0; i < probData[0].size(); i++)
    {
        Scenario s;
        s.probability = probData[0][i];
        s.demandList = vector<Demand>();
        for (int j = 0; j < demandData.size(); j++)
        {
            int sid = Round(demandData[j][2]);
            if (sid != i)
                continue;
            Demand d;
            d.originID = Round(demandData[j][0]);
            d.destinationID = Round(demandData[j][1]);
            d.weight = demandData[j][3];
            d.volumeRatio = demandData[j][4];
            d.revenue = demandData[j][5];
            d.maxTransTime = demandData[j][6];
			int leg_ind = getLegIndex(d.originID, d.destinationID);
			if(leg_ind!=-1)
			{
				double short_time = MAX_NUM;
				for(auto& cf:legList[leg_ind].confList)
				{
					if (short_time > cf.legTime)
						short_time = cf.legTime;
				}
				if (short_time != MAX_NUM && short_time > d.maxTransTime)
					continue;
			}
            s.demandList.push_back(d);
        }
		if (s.demandList.size()>0)
			scenarioList.push_back(s);
    }
	scenarioNum = probData[0].size();

    

    executor = new threadpool(THREAD_POOL_SIZE);
}
