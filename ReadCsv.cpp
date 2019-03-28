#include <fstream>
#include <iostream>
#include <sstream>
#include "ReadCsv.h"

vector<vector<double>> ReadCsv::read(std::string location)
{
	vector <vector <double> > data;
	ifstream infile(location);

	while (infile)
	{
		std::string s;
		if (!getline(infile, s)) break;

		istringstream ss(s);
		vector <double> record;

		while (ss)
		{
			string val;
			if (!getline(ss, val, ',')) break;
			record.push_back(stod(val));
		}

		data.push_back(record);
	}
	return data;
}
