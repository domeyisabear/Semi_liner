#pragma once
#include "SubMIPProblem.h"
#include "InputParam.h"

class SMP_extensive :
	public SubMIPProblem
{
private:
	int sind;
	double lowerBound = -MAX_NUM;

public:
	SMP_extensive(int scenarioIndex);
	SubMIPSolution getOptimalSolution(MasterSolution &ms) override;
	double getLowerBound() override;
	~SMP_extensive()
	{};
};

