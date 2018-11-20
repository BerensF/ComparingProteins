/*

 * second_lower_bound.h

 *

 *  Created on: 26.03.2018

 *      Author: Berens

 *       Function to compute a lower bound of the Gromov Wasserstein Metric for discrete points,
 *          see Section 3.4.2 and 4.2.2 "Third Lower Bound" in "Quantitative comparison of protein isosurfaces 
 *			with approximated Gromov-Wasserstein-distance" from Felix Berens.

 */



#ifndef SECOND_LOWER_BOUND_H_

#define SECOND_LOWER_BOUND_H_



#include "../mathFunctions.h"

#include <vector>

#include <algorithm>



using namespace std;




void list_of_unique_dist(vector< vector<double> > &d1, vector< vector<double> > &d2, vector<double> &listDist)

{

	/* list_of_unique_dist gets two distance matrices an gives a list (listDist) of all occurring distances out. This list is sorted in ascending order */



	for(unsigned int i=0; i < d1.size(); i++) // d1.size() == d2.size(), both of dimension number_of_selected_points x number_of_selected_points

	{

		for(unsigned int j=i;j< d1.size(); j++)

		{

			listDist.push_back(d1[i][j]);

			listDist.push_back(d2[i][j]);

		}

	}



	sort(listDist.begin(), listDist.end());

  	vector<double>::iterator last = unique(listDist.begin(), listDist.end());

	listDist.erase(last, listDist.end());

}



void which_matrix(vector< vector<double> > &d, vector< vector<double> > &list, double t)

{

	/* which_matrix: gets a distance matrix an gives all positions, where the distance is lower or equal as t  */

	
	int count = 0;

	list[0].resize(d.size()*d.size());

	list[1].resize(d.size()*d.size());



	for(unsigned int i = 0; i < d.size(); i++)	//	only the triangular matrix has to be considered, since a distance matrix is symmetric

	{

		for(unsigned int j = i+1; j < d.size(); j++) // without the diagonal

		{

			if(d[i][j] <= t)

			{

				list[0][count] = i;

				list[1][count] = j;

				count++;

				list[0][count] = j;

				list[1][count] = i;

				count++;

			}

		}

		list[0][count] = i;		// the diagonal is always lower or equal t, because on the diagonal only 0 occurs

		list[1][count] = i;

		count++;

	}

	

	list[0].resize(count);

	list[1].resize(count);

}



double distr_of_distances(vector< vector <double> > &d, vector<double> &distribution, double t)

{

	/* distr_of_distances: See "Gromov–Wasserstein Distances and the Metric Approach to Object Matching", Memoli, Definition 5.4 */

	vector< vector<double> >list_index_of_points_distance_le_t(2);

	double muxmu = 0;



	which_matrix(d, list_index_of_points_distance_le_t, t);


	for(unsigned int i = 0; i < list_index_of_points_distance_le_t[0].size(); i++)

	{

		muxmu = muxmu + distribution[list_index_of_points_distance_le_t[0][i]] * distribution[list_index_of_points_distance_le_t[1][i]];

	}



	return(muxmu);

}



double slb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, vector<double> &distribution2, int number_of_selected_points)

{

	/* Here the lower bound is calculated */



	vector < vector<double> > d1(number_of_selected_points, vector<double>(number_of_selected_points,0));

	vector < vector<double> > d2(number_of_selected_points, vector<double>(number_of_selected_points,0));

	vector <double> listDist;

	double t = -1;

	double out = 0;

	const clock_t begin_time0 = clock();

	dist_matrix(d1, points1);

	dist_matrix(d2, points2);



	list_of_unique_dist(d1,d2,listDist);

	

	for(unsigned int i = 0; i < listDist.size()-1; i++)

	{

		t = listDist[i];

		out += (listDist[i+1]-listDist[i])*abs(distr_of_distances(d1, distribution1, t) - distr_of_distances(d2, distribution2, t));

	}



	return(0.5*out);

}

#endif /* SECOND_LOWER_BOUND_H_ */
