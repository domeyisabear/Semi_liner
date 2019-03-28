#include "MP_cutplane_multi.h"
#include "SMP_extensive.h"
#include <iomanip>
#include <fstream>

MP_cutplane_multi::MP_cutplane_multi()
{
	slbp_vector = vector<SubLBProblem*>(param->scenarioNum);
	smp_vector = vector<SubMIPProblem*>(param->scenarioNum);
	std::vector< std::future<void> > results;
	for (int w = 0; w < param->scenarioNum; w++)
	{
		results.push_back(
			param->executor->commit([](vector<SubLBProblem*> *slbp_v, vector<SubMIPProblem*> *smp_v, int i)
		{
			(*slbp_v)[i] = new SubLBProblem(i); (*smp_v)[i] = new SMP_extensive(i);
		}, &slbp_vector, &smp_vector, w)
			);
	}
	for (int w = 0; w < param->scenarioNum;w++)
	{
		results[w].get();
	}
}

MasterSolution MP_cutplane_multi::getSolution()
{
	clock_t start = clock();

	double LB = 0;
	double UB = MAX_NUM;

	MasterSolution re_ms;
	bool violate = false;

	vector<double> subLBs = vector<double>(param->scenarioNum);
	{   // parallelize
		std::vector< std::future<double> > results;
		for (int w = 0; w < param->scenarioNum; w++)
		{
			results.push_back(param->executor->commit([](SubMIPProblem* smp) {return smp->getLowerBound(); }, smp_vector[w]));
		}
		for (int w = 0; w < param->scenarioNum; w++)
		{
			subLBs[w] = results[w].get();
			LB += param->scenarioList[w].probability*subLBs[w];
		}
	}

	double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "UB: " << UB << ", LB:" << LB << ", time elapsed: " << time_elps << endl;
	string workDir = param->getWorkDir();
	ofstream outputfile(workDir + "UB_LB_cutplane_multi_" + std::to_string(dualUsed) + ".csv");

	////////////////// the strengthened method based on lagrangian relaxation ///////////////////
	double rau = 0;// step size for dual price
	double step_size_coeff = 0.8; // coefficient for step size within (0,2)
	vector<vector<double>> omega_list = vector<vector<double>>(param->scenarioNum);	// dual value
	for (int i = 0; i < omega_list.size(); i++)
		omega_list[i] = vector<double>(param->fixPrioList.size());
	vector<vector<int>> a_list = vector<vector<int>>(param->scenarioNum); // solution
	for (int i = 0; i < a_list.size(); i++)
		a_list[i] = vector<int>(param->fixPrioList.size());
	vector<double> a_avg = vector<double>(param->fixPrioList.size()); // solution average
	/////////////////////////////////////////////////////////////////////////////////////////////

	bool first = true;
	int iter_num = 0;

	while (!violate)
	{
		for (int w = 0; w < param->scenarioNum; w++)
		{
			if (dualUsed)
			{
				// update dual value
				for (int i = 0; i < param->fixPrioList.size(); i++)
					omega_list[w][i] += rau*(a_list[w][i] - a_avg[i]);

				// update formulation based on new dual value
				slbp_vector[w]->setDualPenalty(omega_list[w]);
			}

			// calculate the obj value (or LB) and solution for each lower bound subproblem
			std::cout << "problem : " << w << endl;
			vector<MasterSolution> ms_vector = slbp_vector[w]->getSolutions();

			// check feasible
			if (ms_vector.empty()) // infeasible
			{
				std::cout << "UB: " << UB << ", LB:" << MAX_NUM << ", SolN: " << 0 << endl;
				outputfile << setprecision(16) << UB << "," << MAX_NUM << "," << 0 << "\n";
				violate = true;
				break;
			}

			// update lowerbound
			if (subLBs[w] < ms_vector[0].objLB)
			{
				LB = LB - param->scenarioList[w].probability*(subLBs[w] - ms_vector[0].objLB);
				subLBs[w] = ms_vector[0].objLB;
			}

			// check convergence
			if (UB<LB)
			{
				std::cout << "UB: " << UB << ", LB:" << LB << ", SolN: " << 0 << endl;
				outputfile << setprecision(16) << UB << "," << LB << "," << 0 << "\n";
				violate = true;
				break;
			}

			// check time limit
			double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
			if (time_elps>timeLimit)
			{
				std::cout << "UB: " << UB << ", LB:" << LB << ", SolN: " << 0 <<", time elapsed: " << time_elps << endl;
				outputfile << setprecision(16) << UB << "," << LB << "," << 0 << "\n";
				violate = true;
				break;
			}

			// record solution
			vector<int>().swap(a_list[w]);
			a_list[w] = ms_vector[0].a_ij;

			// remove these solutions from all lower bound subproblems
			{
				std::future<void> result;
				for (int i = 0; i < param->scenarioNum; i++)
				{
					result = param->executor->commit([](SubLBProblem* slbp, vector<MasterSolution> &ms_v) {slbp->removeSolutions(ms_v); }, slbp_vector[i], ms_vector);
				}
				result.wait();
			}

			// calculate the obj value (negative profit) of these solutions and update upper bound
			double min_val=MAX_NUM;
			int SolN = ms_vector.size();

			double tmp = 0;
			for (int i = 0; i < param->scenarioNum;i++)
			{
				if (i == w)
					continue;
				tmp += param->scenarioList[i].probability*subLBs[i];
			}
			double level = (UB - tmp) / param->scenarioList[w].probability;

			for (int i = 0; i < ms_vector.size();i++)
			{
				std::vector< std::future<double> > results = std::vector< std::future<double> >(param->scenarioNum);
				for (int j = 0; j < param->scenarioNum; j++)
				{
					results[j] = param->executor->commit([](SubMIPProblem* smp, MasterSolution &ms_) {return smp->getOptimalSolution(ms_).objVal; }, smp_vector[j], ms_vector[i]);
				}

				double objVal = results[w].get();

				if (i > 0)
				{
					objVal = smp_vector[w]->getOptimalSolution(ms_vector[i]).objVal;
					double dualVal = 0;
					for (int j = 0; j < param->fixPrioList.size(); j++)
					{
						dualVal += omega_list[w][j] * ms_vector[i].a_ij[j];
					}
					if (objVal + dualVal >= level) // the obj value of original problem must be larger than UB, so do not need to explore it
					{
						cout << "solution omitted!!" << endl;
						continue;
					}
				}

				double val = param->scenarioList[w].probability*objVal;
				
				for (int j = 0; j < param->scenarioNum; j++)
				{
					if (j == w)
						continue;
					val += param->scenarioList[j].probability*results[j].get();
				}

				if (val<UB)
				{
					UB = val;
					re_ms = ms_vector[i];
					re_ms.objVal = val;
				}
				if (val < min_val)
					min_val = val;
			}

			// check convergence
			if (UB<LB)
			{
				std::cout << "UB: " << UB << ", LB:" << LB << ", SolN: " << SolN << ", time elapsed: " << time_elps << endl;
				outputfile << setprecision(16) << UB << "," << LB << "," << SolN << "\n";
				violate = true;
				break;
			}

			// check time limit
			time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
			if (time_elps>timeLimit)
			{
				std::cout << "UB: " << UB << ", LB:" << LB << ", SolN: " << SolN << ", time elapsed: " << time_elps << endl;
				outputfile << setprecision(16) << UB << "," << LB << "," << SolN << "\n";
				violate = true;
				break;
			}

			std::cout << "UB: " << UB << ", LB:" << LB << ", SolN: " << SolN << ", time elapsed: " << time_elps << endl;
			outputfile << setprecision(16) << UB << "," << LB << "," << SolN << "\n";

			// check iteration limit
			iter_num++;
			if (iterationLim > 0 && iter_num >= iterationLim)
			{
				violate = true;
				break;
			}

			if (dualUsed)
			{
				if (!first || (first && w == param->scenarioNum - 1))
				{
					// calculate average solution
					for (int i = 0; i < param->fixPrioList.size(); i++)
					{
						a_avg[i] = 0;
						for (int j = 0; j < param->scenarioNum; j++)
						{
							a_avg[i] += param->scenarioList[j].probability*a_list[j][i];
						}
					}

					// update step size
					double sum = 0;
					for (int k = 0; k < param->scenarioNum; k++)
					{
						for (int i = 0; i < param->fixPrioList.size(); i++)
						{
							sum += (a_list[k][i] - a_avg[i])*(a_list[k][i] - a_avg[i]);
						}
					}
					rau = step_size_coeff*(min_val - LB) / sqrt(sum);
					std::cout << "step size = " << rau << endl;
				}
			}
		}

		first = false;

		if (violate)
			break;

		outputfile.flush();
	}

	outputfile.close();

	re_ms.objLB = LB<UB ? LB : UB;

	return re_ms;
}
