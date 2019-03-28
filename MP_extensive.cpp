#include "MP_extensive.h"
#include "gurobi_c++.h"
#include "InputParam.h"
#include "tools.h"
#include <fstream>
#include <iomanip>

MP_extensive::MP_extensive()
{
	// initialize variables
	a_ij = addVars(model, 0, 1, GRB_BINARY, "a", param->fixPrioList.size());
	b_w_ij = addVars(model, 0, 1, GRB_BINARY, "b", param->scenarioNum, param->prioList.size());
	x_w_ij_s = GRBVar3d(param->scenarioNum);
	for (int i = 0; i < param->scenarioNum; i++)
	{
		x_w_ij_s[i] = GRBVar2d(param->legNum);
		for (int j = 0; j < param->legNum; j++)
		{
			string str = "x_" + to_string(i) + "_" + to_string(j);
			x_w_ij_s[i][j] = addVars(model, 0, 1, GRB_BINARY, str, param->legList[j].confList.size());
		}
	}
	y_w_i = addVars(model, 0, 1, GRB_BINARY, "y", param->scenarioNum, param->addPortIndices.size());
	g_w_h = GRBVar2d(param->scenarioNum);
	gama_w_h = GRBVar2d(param->scenarioNum);
	for (int i = 0; i < param->scenarioNum; i++)
	{
		string str = "g_" + to_string(i);
		g_w_h[i] = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, str, param->scenarioList[i].demandList.size());
		gama_w_h[i] = addVars(model, 0, 1, GRB_BINARY, "gama_" + to_string(i),param->scenarioList[i].demandList.size());
	}
	f_w_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "f", param->scenarioNum, param->legNum);
	fv_w_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "fv", param->scenarioNum, param->legNum);
	t_w_i = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "t", param->scenarioNum, param->portNum);

	// constraints

	// constraint (1)
	{
		set<set<int>> S;
		for (int i = 0; i < param->fixPortIndices.size(); i++)
		{
			for (int j = 0; j < param->fixPortIndices.size(); j++)
			{
				if (i == j)
					continue;
				int ind_i = param->fixPortIndices[i];
				int ind_j = param->fixPortIndices[j];

				int port1 = param->portList[ind_i].portID;
				int port2 = param->portList[ind_j].portID;
				int index1 = param->getFixPrioIndex(port1, port2);
				int index2 = param->getFixPrioIndex(port2, port1);

				set<int> s;
				s.insert(index1);
				s.insert(index2);

				int len1 = S.size();
				S.insert(s);
				int len2 = S.size();
				if (len2 > len1)
					model.addConstr(a_ij[index1] + a_ij[index2] == 1);
			}
		}
	}

	// constraint (2)
	{
		set<set<int>> S;
		for (int i = 0; i < param->fixPortIndices.size(); i++)
		{
			for (int j = 0; j < param->fixPortIndices.size(); j++)
			{
				if (i == j)
					continue;
				for (int k = 0; k < param->fixPortIndices.size(); k++)
				{
					if ((i == k) || (j == k))
						continue;
					int ind_i = param->fixPortIndices[i];
					int ind_j = param->fixPortIndices[j];
					int ind_k = param->fixPortIndices[k];

					int port1 = param->portList[ind_i].portID;
					int port2 = param->portList[ind_j].portID;
					int port3 = param->portList[ind_k].portID;
					int index1 = param->getFixPrioIndex(port1, port2);
					int index2 = param->getFixPrioIndex(port2, port3);
					int index3 = param->getFixPrioIndex(port3, port1);

					set<int> s;
					s.insert(index1);
					s.insert(index2);
					s.insert(index3);

					int len1 = S.size();
					S.insert(s);
					int len2 = S.size();

					if (len2 > len1)
						model.addConstr(a_ij[index1] + a_ij[index2] + a_ij[index3] <= 2);
				}
			}
		}
	}

	// constraint (3)
	for (int w = 0; w < param->scenarioNum; w++)
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
				model.addConstr(b_w_ij[w][i] + b_w_ij[w][prio2] == 1);
		}
	}

	// constraint (4)
	for (int w = 0; w < param->scenarioNum; w++)
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
						model.addConstr(b_w_ij[w][prio_ij] + b_w_ij[w][prio_jk] + b_w_ij[w][prio_ki] <= 2);
				}
			}
		}
	}

	// constraint (6)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			int id_i = param->fixPrioList[i].first;
			int id_j = param->fixPrioList[i].second;

			int prio2 = param->getPrioIndex(id_i, id_j);
			if (prio2 != -1)
				model.addConstr(b_w_ij[w][prio2] - a_ij[i] == 0);
		}
	}

	// constraint (8)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int a = 0; a < param->legNum; a++)
		{
			GRBLinExpr expr = 0;
			for (GRBVar& var : x_w_ij_s[w][a])
			{
				expr += var;
			}
			model.addConstr(expr <= 1);
		}
	}

	// constraint (9)
	for (int w = 0; w < param->scenarioNum; w++)
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
					for (GRBVar& var : x_w_ij_s[w][in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_w_ij_s[w][out_index])
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
					for (GRBVar& var : x_w_ij_s[w][out_index])
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
					for (GRBVar& var : x_w_ij_s[w][in_index])
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
					for (GRBVar& var : x_w_ij_s[w][in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_w_ij_s[w][out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr - y_w_i[w][i] == 0);
			model.addConstr(in_expr - y_w_i[w][i] == 0);
		}
	}

	// constraint (10)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;
			int leg_ind = param->getLegIndex(id_i, id_j);

			if (leg_ind != -1)
			{
				for (auto& var : x_w_ij_s[w][leg_ind])
				{
					model.addConstr(var - b_w_ij[w][i] <= 0);
				}
			}
		}
	}

	// constraint (11)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		GRBLinExpr expr = 0;
		for (int a = 0; a < param->legNum; a++)
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				int tail_ind = param->getPortIndex(param->legList[a].tailID);
				double voyTime = param->legList[a].confList[s].legTime + param->portList[tail_ind].portTime;
				expr += voyTime * x_w_ij_s[w][a][s];
			}
		model.addConstr(expr <= param->maxVoyageTime);
	}

	// constraint (12) (13)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			GRBLinExpr expr = 0;
			for (auto& var : x_w_ij_s[w][i])
				expr += var;
			model.addConstr(f_w_ij[w][i] - param->shipWeightCapacity * expr <= 0);
			model.addConstr(fv_w_ij[w][i] - param->shipVolumeCapacity * expr <= 0);
		}
	}

	// constraint (18) (19)
	for (int w = 0; w < param->scenarioNum; w++)
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
					wlhs += f_w_ij[w][a];
					vlhs += fv_w_ij[w][a];
				}

				if (tailID == port_i)
				{
					wlhs -= f_w_ij[w][a];
					vlhs -= fv_w_ij[w][a];
				}
			}

			GRBLinExpr wrhs = 0;
			GRBLinExpr vrhs = 0;
			for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
			{
				int originID = param->scenarioList[w].demandList[k].originID;
				int destID = param->scenarioList[w].demandList[k].destinationID;

				if (destID == port_i)
				{
					wrhs += g_w_h[w][k];
					vrhs += param->scenarioList[w].demandList[k].volumeRatio * g_w_h[w][k];
				}

				if (originID == port_i)
				{
					wrhs -= g_w_h[w][k];
					vrhs -= param->scenarioList[w].demandList[k].volumeRatio * g_w_h[w][k];
				}
			}
			model.addConstr(wlhs - wrhs == 0);
			model.addConstr(vlhs - vrhs == 0);
		}
	}

	// constraint (20)
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
		{
			int originID = param->scenarioList[w].demandList[k].originID;
			int destID = param->scenarioList[w].demandList[k].destinationID;

			double demandAmount = param->scenarioList[w].demandList[k].weight;

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
					model.addConstr(g_w_h[w][k] - demandAmount * y_w_i[w][y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_w_h[w][k] - demandAmount <= 0);
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
					model.addConstr(g_w_h[w][k] - demandAmount * y_w_i[w][y_ind] <= 0);
				}
				else
				{
					model.addConstr(g_w_h[w][k] - demandAmount <= 0);
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				model.addConstr(g_w_h[w][k] - demandAmount * b_w_ij[w][prioInd] <= 0);

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
					model.addConstr(g_w_h[w][k] - demandAmount * y_w_i[w][y_origin_ind] <= 0);
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
					model.addConstr(g_w_h[w][k] - demandAmount * y_w_i[w][y_dest_ind] <= 0);
				}
			}
		}
	}

	// constraints (21)~(23)
	for(int w=0;w<param->scenarioNum;w++)
	{
		// constraint (21)
		for(int i=0;i<param->portNum;i++)
		{
			if(i==param->baseStartIndex)
				continue;
			for(int j=0;j<param->portNum;j++)
			{
				if (j == param->baseStartIndex)
					continue;
				if (i == j)
					continue;

				int i_id = param->portList[i].portID;
				int j_id = param->portList[j].portID;

				GRBLinExpr expr1 = 0;
				int leg_ind1 = param->getLegIndex(i_id, j_id);
				if(leg_ind1!=-1)
				{
					for(int s=0;s<param->legList[leg_ind1].confList.size();s++)
					{
						double t_ij = param->legList[leg_ind1].confList[s].legTime + param->portList[i].portTime;
						expr1 += (param->maxVoyageTime + t_ij)*x_w_ij_s[w][leg_ind1][s];
					}
				}

				GRBLinExpr expr2 = 0;
				int leg_ind2 = param->getLegIndex(j_id, i_id);
				if(leg_ind2!=-1)
				{
					for (int s = 0; s < param->legList[leg_ind2].confList.size(); s++)
					{
						double t_ji = param->legList[leg_ind2].confList[s].legTime + param->portList[j].portTime;
						expr2 += (param->maxVoyageTime - t_ji)*x_w_ij_s[w][leg_ind2][s];
					}
				}

				model.addConstr(t_w_i[w][i] - t_w_i[w][j] + expr1 + expr2 <= param->maxVoyageTime);
				//model.addConstr(t_w_i[w][i] - t_w_i[w][j] + expr1 <= param->maxVoyageTime);
			}
		}

		// constraint (22)
		GRBLinExpr expr1 = 0;
		for(int a=0;a<param->legNum;a++)
		{
			int tail_id = param->legList[a].tailID;

			if(tail_id==param->baseStartID)
			{
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					double t_0j = param->legList[a].confList[s].legTime + param->portList[param->baseStartIndex].portTime;
					expr1 += t_0j * x_w_ij_s[w][a][s];
				}
			}
		}
		for(int i=0;i<param->portNum;i++)
		{
			if(i==param->baseStartIndex)
				continue;

			int i_id = param->portList[i].portID;

			GRBLinExpr expr2 = 0;
			for(int a=0;a<param->legNum;a++)
			{
				int tail_id = param->legList[a].tailID;
				int head_id = param->legList[a].headID;
				if(head_id==i_id && tail_id!=param->baseStartID)
				{
					int j_ind = param->getPortIndex(tail_id);
					for(int s=0;s<param->legList[a].confList.size();s++)
					{
						double t_ji = param->legList[a].confList[s].legTime + param->portList[j_ind].portTime;
						expr2 += t_ji * x_w_ij_s[w][a][s];
					}
				}
			}
			model.addConstr(t_w_i[w][i] - expr1 - expr2 >= 0);
			model.addConstr(t_w_i[w][i] - param->maxVoyageTime <= 0);
		}

		// constraint (23)
		model.addConstr(t_w_i[w][param->baseStartIndex] == 0);
	}

	// constraints (24)(25)
	for(int w=0;w<param->scenarioNum;w++)
	{
		for(int h=0;h<param->scenarioList[w].demandList.size();h++)
		{
			double weight = param->scenarioList[w].demandList[h].weight;
			double timeLimit = param->scenarioList[w].demandList[h].maxTransTime;

			int o_ind = param->getPortIndex(param->scenarioList[w].demandList[h].originID);
			int d_ind = param->getPortIndex(param->scenarioList[w].demandList[h].destinationID);

			model.addConstr(param->maxVoyageTime*(gama_w_h[w][h] - 1) - t_w_i[w][d_ind] + t_w_i[w][o_ind] <= 0);
			model.addConstr(t_w_i[w][d_ind] - t_w_i[w][o_ind] - timeLimit - param->maxVoyageTime*(1 - gama_w_h[w][h]) <= 0);
			model.addConstr(g_w_h[w][h] - weight * gama_w_h[w][h] <= 0);
	
			expr_v.push_back(t_w_i[w][d_ind] - t_w_i[w][o_ind]);
		}
	}

	// cuts
	for(int w=0;w<param->scenarioNum;w++)
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
				y[i] += y_w_i[w][ind];
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
					model.addConstr(y[i] - y[j] - b_w_ij[w][prioInd] <= 0);
					model.addConstr(b_w_ij[w][prioInd] - 1 - y[i] + y[j] <= 0);
					if (param->portList[i].portID < param->portList[j].portID)
						model.addConstr(1 - y[i] - y[j] - b_w_ij[w][prioInd] <= 0);
				}

			}
		}
	}

	// objective
	GRBLinExpr obj = 0;
	for (int w = 0; w < param->scenarioNum; w++)
	{
		GRBLinExpr sub_obj = 0;
		for (int i = 0; i < param->legNum; i++)
		{
			int tail_ind = param->getPortIndex(param->legList[i].tailID);
			for (int s = 0; s < param->legList[i].confList.size(); s++)
			{
				double c_ij = param->legList[i].confList[s].legCost + param->portList[tail_ind].portCost;
				sub_obj += c_ij * x_w_ij_s[w][i][s];
			}
		}
		for (int k = 0; k < param->scenarioList[w].demandList.size(); k++)
		{
			sub_obj -= param->scenarioList[w].demandList[k].revenue * g_w_h[w][k];
		}
		obj += param->scenarioList[w].probability * sub_obj;
	}

	model.setObjective(obj, GRB_MINIMIZE);

	isSolved = false;
	model.update();

	int thread_num = std::thread::hardware_concurrency();
	if (thread_num)
		model.set(GRB_IntParam_Threads, thread_num);

	//model.write("model.lp");
	ofstream outputfile(param->getWorkDir() + "problem size.txt");
	outputfile << "binary var : " << model.get(GRB_IntAttr_NumBinVars) << "\n";
	outputfile << "continous var : " << model.get(GRB_IntAttr_NumVars) - model.get(GRB_IntAttr_NumBinVars) << "\n";
	outputfile << "constr : " << model.get(GRB_IntAttr_NumConstrs) << "\n";
	outputfile.close();
}

MasterSolution MP_extensive::getSolution()
{
	if (!isSolved)
	{
		model.set(GRB_DoubleParam_TimeLimit, timeLimit);
		//model.set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();
		isSolved = true;
	}
	//cout << "Optimal value: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
	double objVal = 0;
	if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
	{
		objVal = MAX_NUM;
		ofstream outputfile(param->getWorkDir() + "UB_LB_extensive.csv");
		outputfile << setprecision(16) << objVal << "," << model.get(GRB_DoubleAttr_ObjBound)
			<< "\n";
		outputfile.close();

		MasterSolution ms;
		ms.objVal = MAX_NUM;
		ms.objLB = model.get(GRB_DoubleAttr_ObjBound);
		return ms;
	}

	objVal = model.get(GRB_DoubleAttr_ObjVal);

	try {
		ofstream re_a_ij(param->getWorkDir() + "a_ij.csv");
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			double a = a_ij[i].get(GRB_DoubleAttr_X);
			re_a_ij << param->fixPrioList[i].first << "," << param->fixPrioList[i].second << "," << Round(a) << endl;
		}
		re_a_ij.close();

		ofstream re_b_w_ij(param->getWorkDir() + "b_w_ij.csv");
		for (int i = 0; i < param->prioList.size(); i++)
		{
			re_b_w_ij << param->prioList[i].first << "," << param->prioList[i].second;
			for (int w = 0; w < param->scenarioNum; w++)
			{
				double b = b_w_ij[w][i].get(GRB_DoubleAttr_X);
				re_b_w_ij << "," << Round(b);
			}
			re_b_w_ij << endl;
		}
		re_b_w_ij.close();

		ofstream re_x_w_ij_s(param->getWorkDir() + "x_w_ij_s.csv");
		for (int a = 0; a < param->legNum; a++)
		{
			re_x_w_ij_s << param->legList[a].tailID << "," << param->legList[a].headID;
			for (int w = 0; w < param->scenarioNum; w++)
			{
				re_x_w_ij_s << ",";
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					double x = x_w_ij_s[w][a][s].get(GRB_DoubleAttr_X);
					re_x_w_ij_s << "," << Round(x);
				}

				double f = f_w_ij[w][a].get(GRB_DoubleAttr_X);
				double fv = fv_w_ij[w][a].get(GRB_DoubleAttr_X);
				re_x_w_ij_s << "," << f << "," << fv << ",";
			}
			re_x_w_ij_s << endl;
		}
		re_x_w_ij_s.close();

		ofstream re_y_w_i(param->getWorkDir() + "y_w_i.csv");
		for (int i = 0; i < param->portNum; i++)
		{
			re_y_w_i << param->portList[i].portID;

			int ind = -1;
			for (int j = 0; j < param->addPortIndices.size(); j++)
			{
				if (param->addPortIndices[j] == i)
				{
					ind = j;
					break;
				}
			}
			for (int w = 0; w < param->scenarioNum; w++)
			{
				if (ind == -1)
					re_y_w_i << "," << 1;
				else
				{
					double y = y_w_i[w][ind].get(GRB_DoubleAttr_X);
					re_y_w_i << "," << Round(y);
				}
				re_y_w_i << "," << t_w_i[w][i].get(GRB_DoubleAttr_X) << ",";
			}

			re_y_w_i << endl;
		}
		re_y_w_i.close();

		ofstream re_g_w_h(param->getWorkDir() + "g_w_h.csv");
		int num = 0;
		for (int w = 0; w < param->scenarioNum; w++)
		{
			for (int h = 0; h < param->scenarioList[w].demandList.size(); h++)
			{
				//double gama = gama_w_h[w][h].get(GRB_DoubleAttr_X);

				re_g_w_h << param->scenarioList[w].demandList[h].originID << "," << param->scenarioList[w].demandList[h].
					destinationID << "," << w << "," << g_w_h[w][h].get(GRB_DoubleAttr_X) << endl;// << "," << Round(gama) << "," <<
					//expr_v[num].getValue() << endl;
				num++;
			}
		}
		re_g_w_h.close();

		for (int w = 0; w < param->scenarioNum; w++)
		{
			GRBLinExpr expr = 0;
			for (int a = 0; a < param->legNum; a++)
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					int head_ind = param->getPortIndex(param->legList[a].headID);
					double voyTime = param->legList[a].confList[s].legTime + param->portList[head_ind].portTime;
					expr += voyTime * x_w_ij_s[w][a][s];
				}

			cout << "scenario #" << w << ", voyage time:" << expr.getValue() << endl;
		}
	}
	catch(...){}

	ofstream outputfile(param->getWorkDir() + "UB_LB_extensive.csv");
	outputfile << setprecision(16) << objVal << "," << model.get(GRB_DoubleAttr_ObjBound)
		<< "\n";
	outputfile.close();

	MasterSolution ms;
	ms.objVal = objVal;
	try {
		ms.a_ij = vector<int>(param->fixPrioList.size());
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			double a = a_ij[i].get(GRB_DoubleAttr_X);
			ms.a_ij[i] = Round(a);
		}
	}
	catch(...){}
	return ms;
}
