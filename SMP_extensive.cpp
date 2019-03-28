#include "SMP_extensive.h"
#include "tools.h"
#include <algorithm>

SMP_extensive::SMP_extensive(int scenarioIndex) :sind(scenarioIndex)
{}

SubMIPSolution SMP_extensive::getOptimalSolution(MasterSolution& ms)
{
	set<int> removeLeg = set<int>(); // remove some unnecessary legs according to ms (index)
	set<int> removeDemand = set<int>(); // remove some unnecessary demands according to ms (index)

	// remove unnecessary legs
	vector<int> seq = vector<int>(param->portNum);
	for (int i = 0; i < param->fixPrioList.size(); i++)
	{
		int id_i = param->fixPrioList[i].first;
		seq[param->getPortIndex(id_i)] += ms.a_ij[i];
	}

	{
		for (int i = 0; i < param->fixPortIndices.size(); i++)
		{
			int ind_i = param->fixPortIndices[i];
			for (int j = 0; j < param->fixPortIndices.size(); j++)
			{
				if (i == j)
					continue;
				int ind_j = param->fixPortIndices[j];
				int diff = seq[ind_i] - seq[ind_j];
				if (diff >= 2)
				{
					int leg_ind = param->getLegIndex(param->portList[ind_i].portID, param->portList[ind_j].portID);
					if (leg_ind > -1)
						removeLeg.insert(leg_ind);
					leg_ind = param->getLegIndex(param->portList[ind_j].portID, param->portList[ind_i].portID);
					if (leg_ind > -1)
						removeLeg.insert(leg_ind);
				}
				else if (diff == 1)
				{
					int leg_ind = param->getLegIndex(param->portList[ind_j].portID, param->portList[ind_i].portID);
					if (leg_ind > -1)
						removeLeg.insert(leg_ind);;
				}
			}
		}

		int fixPortNum = param->fixPortIndices.size();
		for (int i = 0; i < fixPortNum; i++)
		{
			int ind_i = param->fixPortIndices[i];
			if (seq[ind_i] == fixPortNum - 1) // the first port in fixed set by ms
			{
				int leg_ind = param->getLegIndex(param->portList[ind_i].portID, param->baseEndID);
				if (leg_ind > -1)
					removeLeg.insert(leg_ind);
			}
			else if (seq[ind_i] == 0) // the last port in the fixed set by ms
			{
				int leg_ind = param->getLegIndex(param->baseStartID, param->portList[ind_i].portID);
				if (leg_ind > -1)
					removeLeg.insert(leg_ind);
			}
			else // the middle ports
			{
				int leg_ind = param->getLegIndex(param->baseStartID, param->portList[ind_i].portID);
				if (leg_ind > -1)
					removeLeg.insert(leg_ind);
				leg_ind = param->getLegIndex(param->portList[ind_i].portID, param->baseEndID);
				if (leg_ind > -1)
					removeLeg.insert(leg_ind);
			}
		}
	}

	// remove unnecessary demands
	for (int i = 0; i < param->scenarioList[sind].demandList.size(); i++)
	{
		int o_id = param->scenarioList[sind].demandList[i].originID;
		int d_id = param->scenarioList[sind].demandList[i].destinationID;

		int o_ind = param->getPortIndex(o_id);
		int d_ind = param->getPortIndex(d_id);

		if (seq[o_ind] < seq[d_ind] &&
			count(param->fixPortIndices.begin(), param->fixPortIndices.end(), o_ind) > 0 &&
			count(param->fixPortIndices.begin(), param->fixPortIndices.end(), d_ind) > 0)
			removeDemand.insert(i);
	}

	// build model

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	int demandNum = param->scenarioList[sind].demandList.size();

	GRBVar1d b_ij = addVars(model, 0, 1, GRB_BINARY, "b", param->prioList.size());
	GRBVar2d x_ij_s = GRBVar2d(param->legNum);
	for (int a = 0; a < param->legNum; a++)
	{
		string str = "x_" + to_string(a);
		x_ij_s[a] = addVars(model, 0, 1, GRB_BINARY, str, param->legList[a].confList.size());
	}
	GRBVar1d y_i = addVars(model, 0, 1, GRB_BINARY, "y", param->addPortIndices.size());
	GRBVar1d g_h = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "g", demandNum);
	GRBVar1d gama_h = addVars(model, 0, 1, GRB_BINARY, "gama", demandNum);
	GRBVar1d f_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "f", param->legNum);
	GRBVar1d fv_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "fv", param->legNum);
	GRBVar1d t_i = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "t", param->portNum);

	// constraint (3)
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
				model.addConstr(b_ij[i] + b_ij[prio2] == 1);
		}
	}

	// constraint (4)
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
						model.addConstr(b_ij[prio_ij] + b_ij[prio_jk] + b_ij[prio_ki] <= 2);
				}
			}
		}
	}

	// constraint (6)
	{
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			int id_i = param->fixPrioList[i].first;
			int id_j = param->fixPrioList[i].second;

			int prio2 = param->getPrioIndex(id_i, id_j);
			if (prio2 != -1)
				model.addConstr(b_ij[prio2] - ms.a_ij[i] == 0);
		}
	}

	// constraint (8)
	{
		for (int a = 0; a < param->legNum; a++)
		{
			if (removeLeg.count(a) > 0)
				continue;
			GRBLinExpr expr = 0;
			for (GRBVar& var : x_ij_s[a])
			{
				expr += var;
			}
			model.addConstr(expr <= 1);
		}
	}

	// constraint (9)
	{
		// compulsory port
		for (int i = 0; i < param->fixPortIndices.size(); i++)
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
				if (in_index != -1 && removeLeg.count(in_index) == 0)
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1 && removeLeg.count(out_index) == 0)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr == 1);
			model.addConstr(in_expr == 1);
		}

		// starting base port
		{
			GRBLinExpr out_expr = 0;
			int port_i = param->portList[param->baseStartIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int out_index = param->getLegIndex(port_i, port_j);
				if (out_index != -1 && removeLeg.count(out_index) == 0)
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
			}
			model.addConstr(out_expr == 1);
		}
		{
			// ending base port
			GRBLinExpr in_expr = 0;
			int port_i = param->portList[param->baseEndIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				if (in_index != -1 && removeLeg.count(in_index) == 0)
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
			}
			model.addConstr(in_expr == 1);
		}

		// optional ports
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
				if (in_index != -1 && removeLeg.count(in_index) == 0)
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1 && removeLeg.count(out_index) == 0)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr - y_i[i] == 0);
			model.addConstr(in_expr - y_i[i] == 0);
		}
	}

	// constraint (10)
	{
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;
			int leg_ind = param->getLegIndex(id_i, id_j);

			if (leg_ind != -1 && removeLeg.count(leg_ind) == 0)
			{
				for (auto& var : x_ij_s[leg_ind])
				{
					model.addConstr(var - b_ij[i] <= 0);
				}
			}
		}
	}

	// constraint (11)
	{
		GRBLinExpr expr = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			if (removeLeg.count(a) > 0)
				continue;
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				int tail_ind = param->getPortIndex(param->legList[a].tailID);
				double voyTime = param->legList[a].confList[s].legTime + param->portList[tail_ind].portTime;
				expr += voyTime * x_ij_s[a][s];
			}
		}

		model.addConstr(expr <= param->maxVoyageTime);
	}

	// constraint (12) (13)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			if (removeLeg.count(i) > 0)
				continue;
			GRBLinExpr expr = 0;
			for (auto& var : x_ij_s[i])
				expr += var;
			model.addConstr(f_ij[i] - param->shipWeightCapacity * expr <= 0);
			model.addConstr(fv_ij[i] - param->shipVolumeCapacity * expr <= 0);
		}
	}

	// constraint (18) (19)
	{
		for (int i = 0; i < param->portNum; i++)
		{
			GRBLinExpr wlhs = 0;
			GRBLinExpr vlhs = 0;
			int port_i = param->portList[i].portID;
			for (int a = 0; a < param->legNum; a++)
			{
				if (removeLeg.count(a) > 0)
					continue;

				int tailID = param->legList[a].tailID;
				int headID = param->legList[a].headID;

				if (headID == port_i)
				{
					wlhs += f_ij[a];
					vlhs += fv_ij[a];
				}

				if (tailID == port_i)
				{
					wlhs -= f_ij[a];
					vlhs -= fv_ij[a];
				}
			}

			GRBLinExpr wrhs = 0;
			GRBLinExpr vrhs = 0;
			for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
			{
				if (removeDemand.count(k) > 0)
					continue;

				int originID = param->scenarioList[sind].demandList[k].originID;
				int destID = param->scenarioList[sind].demandList[k].destinationID;

				if (destID == port_i)
				{
					wrhs += g_h[k];
					vrhs += param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}

				if (originID == port_i)
				{
					wrhs -= g_h[k];
					vrhs -= param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}
			}
			model.addConstr(wlhs - wrhs == 0);
			model.addConstr(vlhs - vrhs == 0);
		}
	}

	// constraint (20)
	{
		for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
		{
			if (removeDemand.count(k) > 0)
				continue;

			int originID = param->scenarioList[sind].demandList[k].originID;
			int destID = param->scenarioList[sind].demandList[k].destinationID;

			double demandAmount = param->scenarioList[sind].demandList[k].weight;

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
					model.addConstr(g_h[k] - demandAmount * y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_h[k] - demandAmount <= 0);
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
					model.addConstr(g_h[k] - demandAmount * y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_h[k] - demandAmount <= 0);
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				model.addConstr(g_h[k] - demandAmount * b_ij[prioInd] <= 0);

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
					model.addConstr(g_h[k] - demandAmount * y_i[y_origin_ind] <= 0);
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
					model.addConstr(g_h[k] - demandAmount * y_i[y_dest_ind] <= 0);
				}
			}
		}
	}

	// constraint (21)
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;
		for (int j = 0; j < param->portNum; j++)
		{
			if (j == param->baseStartIndex)
				continue;
			if (i == j)
				continue;

			int i_id = param->portList[i].portID;
			int j_id = param->portList[j].portID;

			GRBLinExpr expr1 = 0;
			int leg_ind1 = param->getLegIndex(i_id, j_id);
			if (leg_ind1 != -1 && removeLeg.count(leg_ind1) == 0)
			{
				for (int s = 0; s < param->legList[leg_ind1].confList.size(); s++)
				{
					double t_ij = param->legList[leg_ind1].confList[s].legTime + param->portList[i].portTime;
					expr1 += (param->maxVoyageTime + t_ij) * x_ij_s[leg_ind1][s];
				}
			}

			GRBLinExpr expr2 = 0;
			int leg_ind2 = param->getLegIndex(j_id, i_id);
			if (leg_ind2 != -1 && removeLeg.count(leg_ind2) == 0)
			{
				for (int s = 0; s < param->legList[leg_ind2].confList.size(); s++)
				{
					double t_ji = param->legList[leg_ind2].confList[s].legTime + param->portList[j].portTime;
					expr2 += (param->maxVoyageTime - t_ji) * x_ij_s[leg_ind2][s];
				}
			}

			model.addConstr(t_i[i] - t_i[j] + expr1 + expr2 <= param->maxVoyageTime);
		}
	}

	// constraint (22)
	GRBLinExpr expr1 = 0;
	for (int a = 0; a < param->legNum; a++)
	{
		if (removeLeg.count(a) > 0)
			continue;

		int tail_id = param->legList[a].tailID;

		if (tail_id == param->baseStartID)
		{
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				double t_0j = param->legList[a].confList[s].legTime + param->portList[param->baseStartIndex].portTime;
				expr1 += t_0j * x_ij_s[a][s];
			}
		}
	}
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;

		int i_id = param->portList[i].portID;

		GRBLinExpr expr2 = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			if (removeLeg.count(a) > 0)
				continue;

			int tail_id = param->legList[a].tailID;
			int head_id = param->legList[a].headID;
			if (head_id == i_id && tail_id != param->baseStartID)
			{
				int j_ind = param->getPortIndex(tail_id);
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					double t_ji = param->legList[a].confList[s].legTime + param->portList[j_ind].portTime;
					expr2 += t_ji * x_ij_s[a][s];
				}
			}
		}
		model.addConstr(t_i[i] - expr1 - expr2 >= 0);
		model.addConstr(t_i[i] - param->maxVoyageTime <= 0);
	}

	// constraint (23)
	model.addConstr(t_i[param->baseStartIndex] == 0);

	// constraints (24)(25)
	for (int h = 0; h < param->scenarioList[sind].demandList.size(); h++)
	{
		if (removeDemand.count(h) > 0)
			continue;

		double weight = param->scenarioList[sind].demandList[h].weight;
		double timeLimit = param->scenarioList[sind].demandList[h].maxTransTime;

		if (timeLimit > param->maxVoyageTime)
			continue;

		int o_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].originID);
		int d_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].destinationID);

		//model.addConstr(param->maxVoyageTime * (gama_h[h] - 1) - t_i[d_ind] + t_i[o_ind] <= 0);
		model.addConstr(t_i[d_ind] - t_i[o_ind] - timeLimit - param->maxVoyageTime * (1 - gama_h[h]) <= 0);
		model.addConstr(g_h[h] - weight * gama_h[h] <= 0);
	}

	// cuts
	{
		vector<GRBLinExpr> y(param->portNum);
		for (int i = 0; i < param->portNum; i++)
		{
			y[i] = 0;
			int ind = -1;
			for (int j = 0; j < param->addPortIndices.size(); j++)
			{
				if (i == param->addPortIndices[j])
				{
					ind = j;
					break;
				}
			}
			if (ind != -1)
				y[i] += y_i[ind];
			else
				y[i] += 1;
		}
		for (int i = 0; i < param->portNum; i++)
		{
			if (i == param->baseStartIndex || i == param->baseEndIndex)
				continue;
			for (int j = 0; j < param->portNum; j++)
			{
				if (i == j)
					break;
				if (j == param->baseStartIndex || j == param->baseEndIndex)
					continue;

				int prioInd = param->getPrioIndex(param->portList[i].portID, param->portList[j].portID);

				if (prioInd != -1)
				{
					model.addConstr(y[i] - y[j] - b_ij[prioInd] <= 0);
					model.addConstr(b_ij[prioInd] - 1 - y[i] + y[j] <= 0);
					if (param->portList[i].portID < param->portList[j].portID)
						model.addConstr(1 - y[i] - y[j] - b_ij[prioInd] <= 0);
				}

			}
		}
	}

	// objective function
	GRBLinExpr obj = 0;
	GRBLinExpr legCost = 0;
	GRBLinExpr portCost = 0;
	for (int i = 0; i < param->legNum; i++)
	{
		if (removeLeg.count(i) > 0)
			continue;
		int tailID = param->legList[i].tailID;

		for (int s = 0; s < param->legList[i].confList.size(); s++)
		{
			double c_ij_s = param->legList[i].confList[s].legCost;
			legCost += c_ij_s * x_ij_s[i][s];
			portCost += param->portList[param->getPortIndex(tailID)].portCost*x_ij_s[i][s];
		}
	}
	GRBLinExpr frtRev = 0;
	GRBLinExpr frtAmt = 0;
	for (int k = 0; k < demandNum; k++)
	{
		if (removeDemand.count(k) > 0)
			continue;
		frtRev += param->scenarioList[sind].demandList[k].revenue*g_h[k];
		frtAmt += g_h[k];
	}
	obj = legCost + portCost - frtRev;

	model.setObjective(obj, GRB_MINIMIZE);

	model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	model.set(GRB_IntParam_Threads, THREAD_PER_SOLVER);

	////////////////////////////////////////////////
	//model.set(GRB_DoubleParam_MIPGap, 0.005);
	////////////////////////////////////////////////

	model.optimize();

	int s = model.get(GRB_IntAttr_Status);

	SubMIPSolution sms;

	if (s != 2)
	{
		cout << "infeasible" << endl;
		sms.objVal = MAX_NUM;
		return sms;
	}

	sms.b_ij = getOptimalValue(b_ij);
	sms.x_ij_s = getOptimalValue(x_ij_s);
	sms.y_i = getOptimalValue(y_i);
	sms.theta = -MAX_NUM;
	sms.objVal = model.get(GRB_DoubleAttr_ObjVal);
	double tmp = 0;
	sms.frtRev = frtRev.getValue();
	sms.legCost = legCost.getValue();
	sms.portCost = portCost.getValue();
	sms.cargoAmount = frtAmt.getValue();

	return sms;
}

double SMP_extensive::getLowerBound()
{
	if (lowerBound > -MAX_NUM + 10)
		return lowerBound;

	// build model
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	int demandNum = param->scenarioList[sind].demandList.size();

	GRBVar1d b_ij = addVars(model, 0, 1, GRB_BINARY, "b", param->prioList.size());
	GRBVar2d x_ij_s = GRBVar2d(param->legNum);
	for (int a = 0; a < param->legNum; a++)
	{
		string str = "x_" + to_string(a);
		x_ij_s[a] = addVars(model, 0, 1, GRB_CONTINUOUS, str, param->legList[a].confList.size());
	}
	GRBVar1d y_i = addVars(model, 0, 1, GRB_CONTINUOUS, "y", param->addPortIndices.size());
	GRBVar1d g_h = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "g", demandNum);
	GRBVar1d gama_h = addVars(model, 0, 1, GRB_BINARY, "gama", demandNum);
	GRBVar1d f_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "f", param->legNum);
	GRBVar1d fv_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "fv", param->legNum);
	GRBVar1d t_i = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "t", param->portNum);

	// constraint (3)
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
				model.addConstr(b_ij[i] + b_ij[prio2] == 1);
		}
	}

	// constraint (4)
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
						model.addConstr(b_ij[prio_ij] + b_ij[prio_jk] + b_ij[prio_ki] <= 2);
				}
			}
		}
	}

	// constraint (8)
	{
		for (int a = 0; a < param->legNum; a++)
		{
			GRBLinExpr expr = 0;
			for (GRBVar& var : x_ij_s[a])
			{
				expr += var;
			}
			model.addConstr(expr <= 1);
		}
	}

	// constraint (9)
	{
		// compulsory port
		for (int i = 0; i < param->fixPortIndices.size(); i++)
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
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr == 1);
			model.addConstr(in_expr == 1);
		}

		// starting base port
		{
			GRBLinExpr out_expr = 0;
			int port_i = param->portList[param->baseStartIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int out_index = param->getLegIndex(port_i, port_j);
				if (out_index != -1)
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
			}
			model.addConstr(out_expr == 1);
		}
		{
			// ending base port
			GRBLinExpr in_expr = 0;
			int port_i = param->portList[param->baseEndIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				if (in_index != -1)
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
			}
			model.addConstr(in_expr == 1);
		}

		// optional ports
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
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr - y_i[i] == 0);
			model.addConstr(in_expr - y_i[i] == 0);
		}
	}

	// constraint (10)
	{
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;
			int leg_ind = param->getLegIndex(id_i, id_j);

			if (leg_ind != -1)
			{
				for (auto& var : x_ij_s[leg_ind])
				{
					model.addConstr(var - b_ij[i] <= 0);
				}
			}
		}
	}

	// constraint (11)
	{
		GRBLinExpr expr = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				int tail_ind = param->getPortIndex(param->legList[a].tailID);
				double voyTime = param->legList[a].confList[s].legTime + param->portList[tail_ind].portTime;
				expr += voyTime * x_ij_s[a][s];
			}
		}

		model.addConstr(expr <= param->maxVoyageTime);
	}

	// constraint (12) (13)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			GRBLinExpr expr = 0;
			for (auto& var : x_ij_s[i])
				expr += var;
			model.addConstr(f_ij[i] - param->shipWeightCapacity * expr <= 0);
			model.addConstr(fv_ij[i] - param->shipVolumeCapacity * expr <= 0);
		}
	}

	// constraint (18) (19)
	{
		for (int i = 0; i < param->portNum; i++)
		{
			GRBLinExpr wlhs = 0;
			GRBLinExpr vlhs = 0;
			int port_i = param->portList[i].portID;
			for (int a = 0; a < param->legNum; a++)
			{
				int tailID = param->legList[a].tailID;
				int headID = param->legList[a].headID;

				if (headID == port_i)
				{
					wlhs += f_ij[a];
					vlhs += fv_ij[a];
				}

				if (tailID == port_i)
				{
					wlhs -= f_ij[a];
					vlhs -= fv_ij[a];
				}
			}

			GRBLinExpr wrhs = 0;
			GRBLinExpr vrhs = 0;
			for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
			{
				int originID = param->scenarioList[sind].demandList[k].originID;
				int destID = param->scenarioList[sind].demandList[k].destinationID;

				if (destID == port_i)
				{
					wrhs += g_h[k];
					vrhs += param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}

				if (originID == port_i)
				{
					wrhs -= g_h[k];
					vrhs -= param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}
			}
			model.addConstr(wlhs - wrhs == 0);
			model.addConstr(vlhs - vrhs == 0);
		}
	}

	// constraint (20)
	{
		for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
		{
			int originID = param->scenarioList[sind].demandList[k].originID;
			int destID = param->scenarioList[sind].demandList[k].destinationID;

			double demandAmount = param->scenarioList[sind].demandList[k].weight;

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
					model.addConstr(g_h[k] - demandAmount * y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_h[k] - demandAmount <= 0);
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
					model.addConstr(g_h[k] - demandAmount * y_i[y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_h[k] - demandAmount <= 0);
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				model.addConstr(g_h[k] - demandAmount * b_ij[prioInd] <= 0);

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
					model.addConstr(g_h[k] - demandAmount * y_i[y_origin_ind] <= 0);
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
					model.addConstr(g_h[k] - demandAmount * y_i[y_dest_ind] <= 0);
				}
			}
		}
	}

	// constraint (21)
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;
		for (int j = 0; j < param->portNum; j++)
		{
			if (j == param->baseStartIndex)
				continue;
			if (i == j)
				continue;

			int i_id = param->portList[i].portID;
			int j_id = param->portList[j].portID;

			GRBLinExpr expr1 = 0;
			int leg_ind1 = param->getLegIndex(i_id, j_id);
			if (leg_ind1 != -1)
			{
				for (int s = 0; s < param->legList[leg_ind1].confList.size(); s++)
				{
					double t_ij = param->legList[leg_ind1].confList[s].legTime + param->portList[i].portTime;
					expr1 += (param->maxVoyageTime + t_ij) * x_ij_s[leg_ind1][s];
				}
			}

			GRBLinExpr expr2 = 0;
			int leg_ind2 = param->getLegIndex(j_id, i_id);
			if (leg_ind2 != -1)
			{
				for (int s = 0; s < param->legList[leg_ind2].confList.size(); s++)
				{
					double t_ji = param->legList[leg_ind2].confList[s].legTime + param->portList[j].portTime;
					expr2 += (param->maxVoyageTime - t_ji) * x_ij_s[leg_ind2][s];
				}
			}

			model.addConstr(t_i[i] - t_i[j] + expr1 + expr2 <= param->maxVoyageTime);
		}
	}

	// constraint (22)
	GRBLinExpr expr1 = 0;
	for (int a = 0; a < param->legNum; a++)
	{
		int tail_id = param->legList[a].tailID;

		if (tail_id == param->baseStartID)
		{
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				double t_0j = param->legList[a].confList[s].legTime + param->portList[param->baseStartIndex].portTime;
				expr1 += t_0j * x_ij_s[a][s];
			}
		}
	}
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;

		int i_id = param->portList[i].portID;

		GRBLinExpr expr2 = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			int tail_id = param->legList[a].tailID;
			int head_id = param->legList[a].headID;
			if (head_id == i_id && tail_id != param->baseStartID)
			{
				int j_ind = param->getPortIndex(tail_id);
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					double t_ji = param->legList[a].confList[s].legTime + param->portList[j_ind].portTime;
					expr2 += t_ji * x_ij_s[a][s];
				}
			}
		}
		model.addConstr(t_i[i] - expr1 - expr2 >= 0);
		model.addConstr(t_i[i] - param->maxVoyageTime <= 0);
	}

	// constraint (23)
	model.addConstr(t_i[param->baseStartIndex] == 0);

	// constraints (24)(25)
	for (int h = 0; h < param->scenarioList[sind].demandList.size(); h++)
	{
		double weight = param->scenarioList[sind].demandList[h].weight;
		double timeLimit = param->scenarioList[sind].demandList[h].maxTransTime;

		if (timeLimit > param->maxVoyageTime)
			continue;

		int o_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].originID);
		int d_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].destinationID);

		model.addConstr(param->maxVoyageTime * (gama_h[h] - 1) - t_i[d_ind] + t_i[o_ind] <= 0);
		model.addConstr(t_i[d_ind] - t_i[o_ind] - timeLimit - param->maxVoyageTime * (1 - gama_h[h]) <= 0);
		model.addConstr(g_h[h] - weight * gama_h[h] <= 0);
	}

	// cuts
	{
		vector<GRBLinExpr> y(param->portNum);
		for (int i = 0; i < param->portNum; i++)
		{
			y[i] = 0;
			int ind = -1;
			for (int j = 0; j < param->addPortIndices.size(); j++)
			{
				if (i == param->addPortIndices[j])
				{
					ind = j;
					break;
				}
			}
			if (ind != -1)
				y[i] += y_i[ind];
			else
				y[i] += 1;
		}
		for (int i = 0; i < param->portNum; i++)
		{
			if (i == param->baseStartIndex || i == param->baseEndIndex)
				continue;
			for (int j = 0; j < param->portNum; j++)
			{
				if (i == j)
					break;
				if (j == param->baseStartIndex || j == param->baseEndIndex)
					continue;

				int prioInd = param->getPrioIndex(param->portList[i].portID, param->portList[j].portID);

				if (prioInd != -1)
				{
					model.addConstr(y[i] - y[j] - b_ij[prioInd] <= 0);
					model.addConstr(b_ij[prioInd] - 1 - y[i] + y[j] <= 0);
					if (param->portList[i].portID < param->portList[j].portID)
						model.addConstr(1 - y[i] - y[j] - b_ij[prioInd] <= 0);
				}

			}
		}
	}

	// objective function
	GRBLinExpr obj = 0;
	GRBLinExpr legCost = 0;
	GRBLinExpr portCost = 0;
	for (int i = 0; i < param->legNum; i++)
	{
		int tailID = param->legList[i].tailID;

		for (int s = 0; s < param->legList[i].confList.size(); s++)
		{
			double c_ij_s = param->legList[i].confList[s].legCost;
			legCost += c_ij_s * x_ij_s[i][s];
			portCost += param->portList[param->getPortIndex(tailID)].portCost*x_ij_s[i][s];
		}
	}
	GRBLinExpr frtRev = 0;
	for (int k = 0; k < demandNum; k++)
	{
		frtRev += param->scenarioList[sind].demandList[k].revenue*g_h[k];
	}
	obj = legCost + portCost - frtRev;

	model.setObjective(obj, GRB_MINIMIZE);

	model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	model.set(GRB_IntParam_Threads, THREAD_PER_SOLVER);

	model.set(GRB_DoubleParam_TimeLimit, 10);

	model.optimize();

	lowerBound = model.get(GRB_DoubleAttr_ObjBound);

	return lowerBound;
}
