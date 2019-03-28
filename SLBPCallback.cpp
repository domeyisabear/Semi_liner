#include "SLBPCallback.h"
#include "basic.h"

SLBPCallback::SLBPCallback(double timeLimit) : time_lim(timeLimit)
{
	start = clock();
	bestObjVal = MAX_NUM;
}

void SLBPCallback::reset()
{
	start = clock();
	bestObjVal = MAX_NUM;
}

void SLBPCallback::callback()
{
	if (where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
	{
		double objval = getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
		if(objval<bestObjVal-1e-5)
		{
			bestObjVal = objval;
			start = clock();
		}
		else
		{
			clock_t finish = clock();
			if ((double)(finish - start) / CLOCKS_PER_SEC > time_lim)
				abort();
		}
	}
}