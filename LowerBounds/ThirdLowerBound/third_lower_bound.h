/*
 * third_lower_bound.h
 *
 *  Created on: 17.05.2018
 *      Author: Berens
 *       Function to compute a lower bound of the Gromov Wasserstein Metric for discrete points,
 *          see Section 3.4.3 and 4.2.3 "Third Lower Bound" in "Quantitative comparison of protein isosurfaces 
 *			with approximated Gromov-Wasserstein-distance" from Felix Berens.
 */

#ifndef THIRD_LOWER_BOUND_H_
#define THIRD_LOWER_BOUND_H_

#include "CPPCode/mathFunctions.h"
#include <vector>
#include <algorithm>
#include "lp_lib.h"

using namespace std;


void list_of_unique_dist(vector< vector<double> > &d1, vector< vector<double> > &d2, vector<double> &listDist);
void which(vector<double> &d, vector<double> &list, double t, vector<double> &distribution);
double local_shape_distr(vector< vector <double> > &d, vector<double> &distribution, int x, double t);
void obj_func(vector<double> &distribution1, vector<double> &distribution2, vector< vector <double> > &d1, vector< vector <double> > &d2,lprec *lp);
void constrains(vector<double> &distribution1, vector<double> &distribution2,lprec *lp);
double tlb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, vector<double> &distribution2, int number_of_selected_points);

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

void which(vector<double> &d, vector<double> &list, double t, vector<double> &distribution)
{
	/* which: gets a distance vector, that has all distances from a given point, and gives all points, where the distance is lower or equal as t  */

	int count = 0;
	list.resize(d.size());

	for(unsigned int i = 0; i < d.size(); i++)
	{
		if(d[i] <= t)
		{
			list[count] = i;
			count++;
		}
	}
	list.resize(count);
}

double local_shape_distr(vector< vector <double> > &d, vector<double> &distribution, int x, double t)
{
	/* local_shape_distr: See "Gromov-Wasserstein Distances and the Metric Approach to Object Matching", Memoli,  */
	vector<double> list;
	double mu = 0;

	which(d[x], list, t, distribution);

	for(unsigned int i = 0; i < list.size(); i++)
	{
		mu += distribution[list[i]];
	}

	return(mu);
}

void obj_func(vector<double> &distribution1, vector<double> &distribution2, vector< vector <double> > &d1, vector< vector <double> > &d2,lprec *lp)
{
	/* setting the objective function of the optimization program */

	int n = distribution1.size()*distribution2.size();
	double t = 0;
	int l = 1;
	double *row = NULL;
	vector <double> listDist;

	row = (double *) malloc((n+1) * sizeof(*row));
	list_of_unique_dist(d1,d2,listDist);

	for(int i = 0; i <= n; i++)
	{
		row[i] = 0;
	}

	for(int i = 0; i < distribution1.size(); i++)
	{
		for(int j = 0; j < distribution2.size(); j++)
		{
			for(int k = 0; k < listDist.size()-1;k++)
			{
				t = listDist[k];
				row[l] += (listDist[k+1] - listDist[k]) * abs(local_shape_distr(d1, distribution1, i, t) - local_shape_distr(d2, distribution2, j, t));
			}
			row[l] *= 0.5;
			l++;
		}
	}
	set_obj_fn(lp,row);

	free(row);
}

void constrains(vector<double> &distribution1, vector<double> &distribution2,lprec *lp)
{
	/* setting the constrains of the optimization program */
	int n = distribution1.size()*distribution2.size();
	double *row = NULL;
	row = (double *) malloc((n+1) * sizeof(*row));

	set_add_rowmode(lp,TRUE);

	for(int i = 0; i < distribution1.size(); i++)	// adds first half of the constrains
	{
		for(int j = 0; j <= n; j++)	// first every column is filled with 0
		{
			row[j] = 0;
		}

		for(int j = 0; j < distribution2.size(); j++) 
		{
			row[j+i*distribution2.size()+1] = 1; // only some columns are filled with 1, for an explanation see section 3.2.1
		}
		add_constraint(lp,row,EQ,distribution1[i]);
	}

	for(int i = 0; i < distribution2.size(); i++)	// adds second half of the constrains
	{
		for(int j = 0; j <= n; j++)	// first every column is filled with 0
		{
			row[j] = 0;
		}

		for(int j = 0; j < n; j = j + distribution2.size())
		{
			row[j+i+1] = 1; // only some columns are filled with 1, for an explanation see section 3.2.1
		}
		add_constraint(lp,row,EQ,distribution2[i]);
	}

	set_add_rowmode(lp,FALSE);

	free(row);
}

double tlb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, vector<double> &distribution2, int number_of_selected_points)
{
	/* Here the lower bound is calculated */

	vector < vector<double> > d1(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2(number_of_selected_points, vector<double>(number_of_selected_points,0));
	lprec *lp;
	double out;

	const clock_t begin_time0 = clock();
	dist_matrix(d1, points1);
	dist_matrix(d2, points2);

	lp = make_lp(0, number_of_selected_points*number_of_selected_points); // creates the linear model
	set_minim(lp); // set object direction to minimum
	set_verbose(lp, IMPORTANT); // Only important messages are reported, like warnings or errors

	obj_func(distribution1, distribution2, d1, d2, lp); // add the objective function to the model

	constrains(distribution1, distribution2, lp); // add constrains to the model

	solve(lp); // solves the model

	out = get_objective(lp);

	return(out);
}
#endif /* THIRD_LOWER_BOUND_H_ */
