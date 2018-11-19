/*
 * random_points.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *	
 *	This function gets two numbers, the first is the number of points in a file,
 *		the second is the number of points for that the bound should be calculated
 */

#ifndef RANDOM_POINTS_H_
#define RANDOM_POINTS_H_

#include <vector>
#include <cstdlib>
#include <time.h>
#include "mathFunctions.h"

using namespace std;

vector<unsigned int> random_points(unsigned int number_of_all_points, unsigned int number_to_selected)
{
	if(number_of_all_points < number_to_selected)
	{
		cerr << "Can't select more points, then there are. Number of all points: " << number_of_all_points << " Number of selected points: " << number_to_selected;
		exit(-1);
	}

	vector<unsigned int> list;
	vector<unsigned int> out;
	unsigned int random_number;

	for(unsigned int i = 0; i < number_of_all_points; i++)
	{
		list.push_back(i);
	}

	for(unsigned int i = 0; i < number_to_selected; i++)
	{
		random_number = rand() % list.size();
		out.push_back(list[random_number]);
		list.erase(list.begin() + random_number);
	}

	return(out);
}
#endif /* RANDOM_POINTS_H_ */
