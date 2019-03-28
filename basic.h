#pragma once
#include <gurobi_c++.h>
#include<vector>
#include<set>
#include<cmath>

#define Round(x) ((int)(x<0?x-0.5:x+0.5))
#define MAX_NUM 1e10
#define GRBVar3d vector<vector<vector<GRBVar> > >
#define GRBVar2d vector<vector<GRBVar> >
#define GRBVar1d vector<GRBVar>

using namespace std;

extern int timeLimit;
extern bool dualUsed;
extern int solutionCount; // specify the maximum number of solutions obtained in one iteration by cutting plane method
extern int iterationLim;

extern int THREAD_POOL_SIZE;// multi-thread algorithm
extern int THREAD_PER_SOLVER;

struct MasterSolution
{
	vector<int> a_ij; 
	vector<double> Z_w;
	double objVal;
	double objLB;
	void end()
	{
		vector<int>().swap(a_ij);
		vector<double>().swap(Z_w);
	}
};

struct MasterSolution_continuous
{
	vector<double> a_ij;
	vector<double> Z_w;
	void end()
	{
		vector<double>().swap(a_ij);
		vector<double>().swap(Z_w);
	}
};

struct SubDualSolution
{
	vector<double> xi_ij;
	double other_const;
	double objVal;
	void end()
	{
		vector<double>().swap(xi_ij);
	}
};

struct SubMIPSolution
{
	vector<vector<int>> x_ij_s;
	vector<int> b_ij;
	vector<int> y_i;
	double theta;
	double objVal;
	double frtRev;
	double portCost;
	double legCost;
	double cargoAmount;

	void end()
	{
		for (auto& tmp : x_ij_s)
			vector<int>().swap(tmp);
		vector<vector<int>>().swap(x_ij_s);
		vector<int>().swap(b_ij);
		vector<int>().swap(y_i);
	}
};

struct SubMIPSolution_continuous
{
	vector<vector<double>> x_ij_s;
	vector<double> b_ij;
	vector<double> y_i;
	double theta;
	double objVal;

	void end()
	{
		for (auto& tmp : x_ij_s)
			vector<double>().swap(tmp);
		vector<vector<double>>().swap(x_ij_s);
		vector<double>().swap(b_ij);
		vector<double>().swap(y_i);
	}
};

struct SubSubDualSolution
{
	vector<double> rau_ij_d;
	double other_const;
	double objVal;
	void end()
	{
		vector<double>().swap(rau_ij_d);
	}
};

struct SubSolution
{
	vector<int> a_ij;
	double objVal;
	vector<SubDualSolution> subDualSolutions;
};

////////////////////////////////
struct SubInCumbentSolution
{
	vector<int> x_ij;
	void end()
	{
		vector<int>().swap(x_ij);
	}
};

class SortIncumSolution
{
public:
	bool operator() (const SubInCumbentSolution &a, const SubInCumbentSolution &b) const
	{
		int len1 = a.x_ij.size();
		int len2 = b.x_ij.size();
		int len = len1 > len2 ? len2 : len1;

		for (int i = 0; i < len; i++)
		{
			if (a.x_ij[i] < b.x_ij[i])
				return true;
			else if (a.x_ij[i] > b.x_ij[i])
				return false;
		}
		if (len2 > len)
			return true;
		return false;
	}
};

