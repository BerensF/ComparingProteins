/*
 * first_lower_bound.h
 *
 *  Created on: 17.10.2018
 *      Author: Berens
 *       Function to compute a lower bound of the Gromov Wasserstein Metric for discrete points,
 *          see Section 3.4.1 and 4.2.1 "First Lower Bound" in "Quantitative comparison of protein isosurfaces 
 *			with approximated Gromov-Wasserstein-distance" from Felix Berens.
 */

#ifndef FIRST_LOWER_BOUND_H_
#define FIRST_LOWER_BOUND_H_

#include "../mathFunctions.h"
#include <vector>
#include <algorithm>

using namespace std;


void list_of_unique_eccent(vector<double> &ecc1, vector<double> &ecc2, vector<double> &listEcc);
double eccentricities(vector< vector <double> > &d, vector<double> &distribution, int i);
void all_eccentricities(vector< vector <double> > &d, vector<double> &distribution,vector<double> &ecc);
void which_index(vector<double> &ecc, vector<double> &list, double u);
double distr_of_eccentricities(vector<double> &ecc, vector<double> &distribution, double u);
double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, vector<double> &distribution2, int number_of_selected_points);

void list_of_unique_eccent(vector<double> &ecc1, vector<double> &ecc2, vector<double> &listEcc)
{
	/* list_of_unique_eccent gets two eccentricities vectors and gives a vector (listDist) of all occurring eccentricities out. This list is sorted in ascending order */
	
	listEcc.resize(ecc1.size()+ecc2.size());
	
	for(unsigned int i=0; i < ecc1.size(); i++)
	{
		listEcc[i] = ecc1[i];
		listEcc[i+ecc1.size()] = ecc2[i];
	}

	sort(listEcc.begin(), listEcc.end());
  	vector<double>::iterator last = std::unique(listEcc.begin(), listEcc.end());
	listEcc.erase(last, listEcc.end());
}

double eccentricities(vector< vector <double> > &d, vector<double> &distribution, int i)
{
	/* eccentricities: See "2.4.1 First Lower Bound" in "Quantitative comparison of protein isosurfaces 
 *			with approximated Gromov-Wasserstein-distance" Definition 2.17 */
	double s = 0;

	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		s += distribution[k]*d[i][k];
	}
	return(s);
}

void all_eccentricities(vector< vector <double> > &d, vector<double> &distribution,vector<double> &ecc)
{
	/* calculates all eccentricities for a given set of points */
	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		ecc[k] = eccentricities(d,distribution,k);
	}
}

void which_index(vector<double> &ecc, vector<double> &list, double u)
{
	// Stores in list all indices of points that have an eccentricity lower than u
	int count = 0;
	list.resize(ecc.size());

	for(unsigned int i = 0; i < ecc.size(); i++)
	{
		if(ecc[i] <= u)
		{
			list[count] = i;
			count++;
		}
	}
	list.resize(count);
}

double distr_of_eccentricities(vector<double> &ecc, vector<double> &distribution, double u)
{
	/* distr_of_eccentricities: Is the probability if one choose a point randomly, with the given distribution, that the eccentricity is lower than u*/
	vector<double>list_index_of_points_ecc_le_u;
	double mu = 0;

	which_index(ecc, list_index_of_points_ecc_le_u, u);

	for(unsigned int i = 0; i < list_index_of_points_ecc_le_u.size(); i++)
	{
		mu = mu + distribution[list_index_of_points_ecc_le_u[i]];
	}
	return(mu);
}


double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, vector<double> &distribution2, int number_of_selected_points)
{
	/* This function calculates the lower bound */
	vector < vector<double> > d1(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector <double>eccentricities1(number_of_selected_points);
	vector <double>eccentricities2(number_of_selected_points);
	vector<double> ecc1(number_of_selected_points);
	vector<double> ecc2(number_of_selected_points);
	vector<double> listEcc(number_of_selected_points+number_of_selected_points);
	double u;
	double out;
	dist_matrix(d1, points1); // calculation of the distance matrices
	dist_matrix(d2, points2);

	all_eccentricities(d1, distribution1,ecc1); // calculation of all eccentricities
	all_eccentricities(d2, distribution2,ecc2);

	list_of_unique_eccent(ecc1, ecc2, listEcc); // mergin and sorting the list of eccentricities

	// calculating the lower bound
	for(unsigned int i = 0; i < listEcc.size()-1; i++)
	{
		u = listEcc[i];
		out += (listEcc[i+1]-listEcc[i])*abs(distr_of_eccentricities(ecc1, distribution1, u) - distr_of_eccentricities(ecc2, distribution2, u));
	}
	return(0.5*out);
}
#endif /* FIRST_LOWER_BOUND_H_ */
