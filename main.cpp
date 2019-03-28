#include "basic.h"
#include <iomanip>
#include <fstream>
#include "MP_cutplane_multi.h"
#include "MP_extensive.h"
#include "MP_cutplane.h"
#include "MP_benders.h"
#include "CJsonObject.hpp"
#include "SMP_extensive.h"

using namespace std;

void solve_cutplane(bool dualStrengthen, double timeLim)
{
	dualUsed = dualStrengthen;
	timeLimit = timeLim;

	clock_t start, finish;
	start = clock();

	MP_cutplane* mp = new MP_cutplane();
	mp->getSolution();

	finish = clock();

	ofstream outputfile(param->getWorkDir() + "cutplane_" + to_string(dualStrengthen) + ".txt");
	outputfile << "solution time : " << setprecision(8) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	outputfile.close();

	cout << "finish!" << endl;

	delete mp;
}

void solve_cutplane_multi(bool dualStrengthen, double timeLim)
{
	dualUsed = dualStrengthen;
	timeLimit = timeLim;

	clock_t start, finish;
	start = clock();

	MP_cutplane_multi* mp = new MP_cutplane_multi();
	mp->getSolution();

	finish = clock();

	ofstream outputfile(param->getWorkDir() + "cutplane_multi_" + to_string(dualStrengthen) + ".txt");
	outputfile << "solution time : " << setprecision(8) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	outputfile.close();

	cout << "finish!" << endl;

	delete mp;
}

void solve_benders(double timeLim)
{
	timeLimit = timeLim;

	clock_t start, finish;
	start = clock();

	MP_benders* mp = new MP_benders();
	mp->getSolution();

	finish = clock();

	ofstream outputfile(param->getWorkDir() + "benders.txt");
	outputfile << "solution time : " << setprecision(8) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	outputfile.close();

	cout << "finish!" << endl;

	delete mp;
}

void solve_extensive(double timeLim)
{
	timeLimit = timeLim;

	clock_t start, finish;
	start = clock();

	MP_extensive* mp = new MP_extensive();
	mp->getSolution();

	finish = clock();

	ofstream outputfile(param->getWorkDir() + "extensive.txt");
	outputfile << "solution time : " << setprecision(8) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	outputfile.close();

	cout << "finish!" << endl;

	delete mp;
}

void save_solution(MasterSolution& ms, vector<SubMIPSolution>& sub_solutions)
{ 
	vector<double> avg_speed = vector<double>(param->scenarioNum);
	ofstream outputroutes(param->getWorkDir() + "result_routes.csv");
	ofstream outputspeeds(param->getWorkDir() + "result_routes_speeds.csv");
	for (int i = 0; i < param->scenarioNum; i++)
	{
		vector<int> route = vector<int>();
		vector<double> route_speed = vector<double>();

		int p1 = param->baseStartID;
		route.push_back(p1);
		
		while (1)
		{
			for (int j = 0; j < param->legNum; j++)
			{
				int tailID = param->legList[j].tailID;
				int headID = param->legList[j].headID;
				bool flag = false;
				if(tailID==p1)
				{
					for (int k = 0; k < param->legList[j].confList.size(); k++)
					{
						if(sub_solutions[i].x_ij_s[j][k] == 1)
						{
							p1 = headID;
							flag = true;
							route_speed.push_back(k + 10);/////////////////////////
							break;
						}
					}
					if (flag)
						break;
				}
			}
			route.push_back(p1);
			if (p1 == param->baseEndID)
				break;
		}
		for (int j = 0; j < route.size() - 1; j++)
		{
			outputroutes << route[j] << ",";
		}
		outputroutes << route[route.size() - 1] << "\n";

		for (int j = 0; j < route_speed.size() - 1; j++)
		{
			outputspeeds << route_speed[j] << ",";
		}
		outputspeeds << route_speed[route_speed.size() - 1] << "\n";

		double tmp_sum = 0;
		for (auto j : route_speed)
			tmp_sum += j;
		avg_speed[i] = tmp_sum / route_speed.size();
	}
	outputroutes.close();
	outputspeeds.close();
	 
	ofstream outputresult(param->getWorkDir() + "result_summary.csv");
	outputresult << "port call cost,ship operation cost,freight revenue,total cargo weight,average speed,profit\n";
	for (int i = 0; i < param->scenarioNum; i++)
	{
		outputresult << setprecision(16) << sub_solutions[i].portCost << "," << sub_solutions[i].legCost << "," <<
			sub_solutions[i].frtRev << "," << sub_solutions[i].cargoAmount << "," << avg_speed[i] << "," << -sub_solutions[i].objVal << "\n";
	}
	outputresult.close();
}

void exe_application_case(bool primalDual, bool multiCut, int timeLim)
{
	dualUsed = primalDual;
	timeLimit = timeLim;

	clock_t start, finish;
	start = clock();

	MasterSolution ms;

	if (!multiCut)
	{
		MP_cutplane* mp = new MP_cutplane();
		ms = mp->getSolution();
		delete mp;
	}
	else
	{
		MP_cutplane_multi* mpm = new MP_cutplane_multi();
		ms = mpm->getSolution();
		delete mpm;
	}

	finish = clock();

	ofstream outputfile(param->getWorkDir() + "solution time " + to_string(multiCut) + ".txt");
	outputfile << "solution time : " << setprecision(8) << (double)(finish - start) / CLOCKS_PER_SEC << "s" << endl;
	outputfile.close();

	vector<SubMIPSolution> sub_solutions = vector<SubMIPSolution>();

	std::vector< std::future<SubMIPSolution> > results;
	for (int i = 0; i < param->scenarioNum; i++)
	{
		results.push_back(param->executor->commit(
			[](MasterSolution& ms_,int i_)
			{
			SMP_extensive* smp = new SMP_extensive(i_);
			SubMIPSolution sms = smp->getOptimalSolution(ms_);
			delete smp;
			return sms;
			},
			ms,i));
	}
	for (int i = 0; i < param->scenarioNum; i++)
		sub_solutions.push_back(results[i].get());

	save_solution(ms, sub_solutions);

	cout << "finish!" << endl;
}

int main(int argc, char *argv[])
{
	string file_name = "conf.json";
	if(argc>1)
	{
		file_name = argv[1];
	}
	
	ifstream json_file(file_name);
	if (json_file.is_open())
	{
		string line;
		string json_str;
		while (getline(json_file, line))
		{
			json_str = json_str + line;
		}

		json_file.close();

		//std::thread::hardware_concurrency();

		neb::CJsonObject oJson(json_str);

		oJson.Get("threadPoolSize", THREAD_POOL_SIZE);
		oJson.Get("threadNumPerSolver", THREAD_PER_SOLVER);
		oJson.Get("solutionCount", solutionCount);
		oJson.Get("timeLimit", timeLimit);
		oJson.Get("iterationLimit", iterationLim);

		cout << "----------------------------------------" << endl;
		cout << "thread pool size : " << THREAD_POOL_SIZE << endl;
		cout << "thread number per solver : " << THREAD_PER_SOLVER << endl;
		cout << "time limit : " << timeLimit << endl;
		cout << "iteration limit : " << iterationLim << endl;
		cout << "solution pool size (multi-cut method) : " << solutionCount << endl;
		cout << "----------------------------------------" << endl;

		int size = oJson["data"].GetArraySize();

		for (int i = 0; i < size; i++)
		{
			string path = "";
			oJson["data"][i].Get("path", path);

			cout << "path:" << path << endl;

			const int tmp_threadpoolsize = THREAD_POOL_SIZE;
			const int tmp_threadnum = THREAD_PER_SOLVER;

			int threadpoolsize = 0;
			oJson["data"][i].Get("threadPoolSize", threadpoolsize);

			int threadnum = 0;
			oJson["data"][i].Get("threadNumPerSolver", threadnum);

			if (threadpoolsize > 0 && threadnum > 0)
			{
				THREAD_POOL_SIZE = threadpoolsize;
				THREAD_PER_SOLVER = threadnum;
			}

			cout << "pool size : " << THREAD_POOL_SIZE << endl;
			cout << "thread per solver : " << THREAD_PER_SOLVER << endl;

			try
			{
				param = new InputParam(path);
			}
			catch (...)
			{
				cout << "problems occur when reading input parameters from path:" << path << endl;
				continue;
			}

			try
			{
				int method = 0;
				oJson["data"][i].Get("calSolution", method);
				if(method)
				{
					int primaldual = 0;
					int multiCut = 0;
					oJson["data"][i].Get("multiCutForAppCase", multiCut);
					oJson["data"][i].Get("dualForAppCase", primaldual);
					exe_application_case(primaldual,multiCut, timeLimit);
				}
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("cutPlane0", method);
				if (method)
					solve_cutplane(false, timeLimit);
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("cutPlane1", method);
				if (method)
					solve_cutplane(true, timeLimit);
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("cutPlaneMulti0", method);
				if (method)
					solve_cutplane_multi(false, timeLimit);
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("cutPlaneMulti1", method);
				if (method)
					solve_cutplane_multi(true, timeLimit);
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("benders", method);
				if (method)
					solve_benders(timeLimit);
			}
			catch (...) {}

			try {
				int method = 0;
				oJson["data"][i].Get("extensive", method);
				if (method)
					solve_extensive(timeLimit);
			}
			catch (...) {}

			delete param;

			THREAD_POOL_SIZE = tmp_threadpoolsize;
			THREAD_PER_SOLVER = tmp_threadnum;
		}
		oJson.Clear();
	}
	else
	{
		cout << "Please provide the configuration file (.json) in arguments." << endl;
	}

	//getchar();
	return 0;
}
