#include "LinerProblem.h"
#include <fstream>
#include <iomanip>

LinerProblem::LinerProblem()
{
	a_ij = addVars(model, 0, 1, GRB_BINARY, "a", param->prioList.size());
	x_ij = addVars(model, 0, 1, GRB_BINARY, "x", param->legNum);
	y_i = addVars(model, 0, 1, GRB_BINARY, "y", param->addPortIndices.size());
	g_w_od = new GRBVar*[param->scenarioNum];
	for (int i = 0; i < param->scenarioNum; i++)
	{
		//char str[10];
        string str="g_"+to_string(i);
		g_w_od[i] = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, str, param->scenarioList[i].demandList.size());
	}
	f_w_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "f", param->scenarioNum, param->legNum);

	// constraints
	// constraint (1)
	{
		set<set<int>> S;
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;

			int prio2 = param->getPrioIndex(id_j, id_i);

			set<int> s;
			s.insert(i);
			s.insert(prio2);

			int len1 = S.size();
			S.insert(s);
			int len2 = S.size();

			if (len2 > len1)
				model.addConstr(a_ij[i] + a_ij[prio2] == 1);
		}
	}

	// constraint (2)
	{
		set<set<int>> S;
		for (int i = 0; i < param->portNum; i++)
		{
			if (i == param->baseStartIndex || i == param->baseEndIndex)
				continue;
			for (int j = 0; j < param->portNum; j++)
			{
				if (j == param->baseStartIndex || j == param->baseEndIndex)
					continue;
				if (i == j)
					continue;
				for (int k = 0; k < param->portNum; k++)
				{
					if (k == param->baseStartIndex || k == param->baseEndIndex)
						continue;
					if (j == k || k == i)
						continue;
					int port_i = param->portList[i].portID;
					int port_j = param->portList[j].portID;
					int port_k = param->portList[k].portID;
					int prio_ij = param->getPrioIndex(port_i, port_j);
					int prio_jk = param->getPrioIndex(port_j, port_k);
					int prio_ki = param->getPrioIndex(port_k, port_i);
					if (prio_ij == -1 || prio_jk == -1 || prio_ki == -1)
						continue;

					set<int> s;
					s.insert(prio_ij);
					s.insert(prio_jk);
					s.insert(prio_ki);

					int len1 = S.size();
					S.insert(s);
					int len2 = S.size();

					if (len2 > len1)
						model.addConstr(a_ij[prio_ij] + a_ij[prio_jk] + a_ij[prio_ki] <= 2);
				}
			}
		}
	}

	// constraints (3) ~ (4)
	{
		for (int i = 0; i < param->fixPortIndices.size(); i++) // compulsory port
		{
			GRBLinExpr in_expr = 0;
			GRBLinExpr out_expr = 0;
			int ind_i = param->fixPortIndices[i];
			int port_i = param->portList[ind_i].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				int out_index = param->getLegIndex(port_i, port_j);
				if (in_index != -1)
					in_expr += x_ij[in_index];
				if (out_index != -1)
					out_expr += x_ij[out_index];
			}
			model.addConstr(out_expr == 1);
			model.addConstr(in_expr == 1);
		}
		GRBLinExpr out_expr = 0;// starting base port
		int port_i = param->portList[param->baseStartIndex].portID;
		for (int j = 0; j < param->portNum; j++)
		{
			int port_j = param->portList[j].portID;
			if (port_i == port_j)
				continue;
			int out_index = param->getLegIndex(port_i, port_j);
			if (out_index != -1)
				out_expr += x_ij[out_index];
		}
		model.addConstr(out_expr == 1);
		GRBLinExpr in_expr = 0; // ending base port
		port_i = param->portList[param->baseEndIndex].portID;
		for (int j = 0; j < param->portNum; j++)
		{
			int port_j = param->portList[j].portID;
			if (port_i == port_j)
				continue;
			int in_index = param->getLegIndex(port_j, port_i);
			if (in_index != -1)
				in_expr += x_ij[in_index];
		}
		model.addConstr(in_expr == 1);
	}

	// constraints (5) ~ (6)
	for (int i = 0; i < param->addPortIndices.size(); i++)
	{
		GRBLinExpr in_expr = 0;
		GRBLinExpr out_expr = 0;
		int port_i = param->portList[param->addPortIndices[i]].portID;
		for (int j = 0; j < param->portNum; j++)
		{
			int port_j = param->portList[j].portID;
			if (port_i == port_j)
				continue;
			int in_index = param->getLegIndex(port_j, port_i);
			int out_index = param->getLegIndex(port_i, port_j);
			if (in_index != -1)
				in_expr += x_ij[in_index];
			if (out_index != -1)
				out_expr += x_ij[out_index];
		}
		model.addConstr(out_expr - y_i[i] == 0);
		model.addConstr(in_expr - y_i[i] == 0);
	}

	// constraints (7)
	for (int i = 0; i < param->prioList.size(); i++)
	{
		int id_i = param->prioList[i].first;
		int id_j = param->prioList[i].second;
		int leg_ind = param->getLegIndex(id_i, id_j);
		if (leg_ind != -1)
			model.addConstr(a_ij[i] - x_ij[leg_ind] >= 0);
	}

	// constraint (11)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			model.addConstr(f_w_ij[w][i] - param->shipCapacity*x_ij[i] <= 0);
		}
	}

	// constraint (12)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->portNum; i++)
		{
			GRBLinExpr lhs = 0;
			int port_i = param->portList[i].portID;
			for (int a = 0; a < param->legNum; a++)
			{
				int tailID = param->legList[a].tailID;
				int headID = param->legList[a].headID;

				if (headID == port_i)
					lhs += f_w_ij[w][a];
				if (tailID == port_i)
					lhs -= f_w_ij[w][a];
			}
			GRBLinExpr rhs = 0;
			for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
			{
				int originID = param->scenarioList[w].demandList[k].originID;
				int destID = param->scenarioList[w].demandList[k].destinationID;

				if (destID == port_i)
					rhs += g_w_od[w][k];

				if (originID == port_i)
					rhs -= g_w_od[w][k];
			}
			model.addConstr(lhs - rhs == 0);
		}
	}
	// constraint (13)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
		{
			int originID = param->scenarioList[w].demandList[k].originID;
			int destID = param->scenarioList[w].demandList[k].destinationID;

			double demandAmount = param->scenarioList[w].demandList[k].amount;

			if (originID == param->baseStartID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(destID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					model.addConstr(g_w_od[w][k] - demandAmount*y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_w_od[w][k] - demandAmount <= 0);
				}
			}
			else if (destID == param->baseEndID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(originID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					model.addConstr(g_w_od[w][k] - demandAmount*y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_w_od[w][k] - demandAmount <= 0);
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				model.addConstr(g_w_od[w][k] - demandAmount*a_ij[prioInd] <= 0);

				int originInd = param->getPortIndex(originID);
				int destInd = param->getPortIndex(destID);

				int y_origin_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == originInd)
					{
						y_origin_ind = i;
						break;
					}
				}
				if (y_origin_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == originInd)
						{
							y_origin_ind = i;
							break;
						}
					}
					model.addConstr(g_w_od[w][k] - demandAmount*y_i[y_origin_ind] <= 0);
				}

				int y_dest_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == destInd)
					{
						y_dest_ind = i;
						break;
					}
				}
				if (y_dest_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == destInd)
						{
							y_dest_ind = i;
							break;
						}
					}
					model.addConstr(g_w_od[w][k] - demandAmount*y_i[y_dest_ind] <= 0);
				}
			}
		}
	}

	// constraint (14)
	{
		GRBLinExpr tmp_expr = 0;
		for (int i = 0; i < param->addPortIndices.size(); i++)
		{
			tmp_expr += y_i[i];
		}
		model.addConstr(tmp_expr <= param->maxAddxPortCallNum);
	}

	// objective
	GRBLinExpr obj = 0;
	for (int i = 0; i < param->legNum; i++)
	{
		int tailID = param->legList[i].tailID;
		int headID = param->legList[i].headID;
		double c_ij = param->legList[i].legCost;
		legCost += c_ij*x_ij[i];
		portCost += (param->portList[param->getPortIndex(tailID)].portCost + param->portList[param->getPortIndex(headID)].portCost) / 2.0 *x_ij[i];
	}
	obj = legCost + portCost;

	for (int w = 0; w < param->scenarioNum; w++)
	{
		GRBLinExpr sub_rev = 0;
		GRBLinExpr sub_amt = 0;
		for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
		{
			sub_rev += param->scenarioList[w].demandList[k].rate*g_w_od[w][k];
			sub_amt += g_w_od[w][k];
		}
		frtRev.push_back(sub_rev);
		frtAmt.push_back(sub_amt);
		obj -= param->scenarioList[w].probability*sub_rev;
	}

	model.setObjective(obj, GRB_MINIMIZE);
}

void LinerProblem::getSolution()
{
	model.optimize();
	vector<int> val_x_ij = getOptimalValue(x_ij, param->legNum);

	//vector<int> val_y_i = getOptimalValue(y_i, param->addPortIndices.size());
	//for (int i = 0; i < val_y_i.size();i++)
	//{
	//	cout << val_y_i[i] << endl;
	//}


	cout << "save results..." << endl;

	string workDir = param->getWorkDir();

	ofstream outputroutes(workDir + "result_routes_liner.csv");
	vector<int> route = vector<int>();
	int p1 = param->baseStartID;
	route.push_back(p1);
	while (1)
	{
		for (int j = 0; j < param->legNum; j++)
		{
			int tailID = param->legList[j].tailID;
			int headID = param->legList[j].headID;
			if (tailID == p1 && val_x_ij[j] == 1)
			{
				p1 = headID;
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
	outputroutes.close();

	ofstream outputresult(workDir + "result_summary_liner.csv");
	outputresult << "port call cost, ship operation cost, freight revenue, total cargo weight, profit\n";
	double leg_cost = legCost.getValue();
	double port_cost = portCost.getValue();
	for (int i = 0; i < param->scenarioNum; i++)
	{
		double rev = frtRev[i].getValue();
		double amt = frtAmt[i].getValue();
		outputresult << std::setprecision(16) << port_cost << "," << leg_cost << "," << rev << "," << amt << "," << rev - port_cost - leg_cost << "\n";
	}
	outputresult.close();

	cout << "finish!" << endl;
}