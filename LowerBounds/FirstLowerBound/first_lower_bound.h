/*
 * first_lower_bound.h
 *
 *  Created on: 17.10.2018
 *      Author: Berens
 *       Function to compute a lower bound of the Gromov Wasserstein Metric for discrete points,
 *          see section 3.4.1 and 4.2.1 "First Lower Bound" in "Quantitative comparison of protein isosurfaces 
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
double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, 
	   vector<double> &distribution2, int number_of_selected_points);

void list_of_unique_eccent(vector<double> &ecc1, vector<double> &ecc2, vector<double> &listEcc)
{
	/*
		list_of_unique_eccent: 
		Gets two eccentricities vectors and gives a vector (listEcc) of all occurring eccentricities out.
		This list is sorted in ascending order.

		Input 
		-----
		ecc1 - vector containing the eccentricities for points of the first object
		ecc2 - vector containing the eccentricities for points of the second object
		listEcc - here all eccentricities occurring eccentricities will be stored

		Output
		------
		No output

		Complexity
    		----------
		O(N·log(N)), where N = std::distance(listEcc.begin(), listEcc.end())

		Example
		-------
		vector<double> listEcc(12);
		vector <double > ecc1{ 5, 10, 13, 2, 2, 2};
		vector <double > ecc2{ 3, 34, 1, 2 ,5, 1};
		
		list_of_unique_eccent(ecc1, ecc2, listEcc);
		
		=> listEcc = {1, 2, 3, 5, 10, 13, 34}
	*/
	listEcc.resize(ecc1.size() + ecc2.size());
	
	for(unsigned int i = 0; i < ecc1.size(); i++)
	{
		listEcc[i] = ecc1[i];
	}
	
	for(unsigned int i = 0; i < ecc2.size(); i++)
	{
		listEcc[i+ecc1.size()] = ecc2[i];
	}

	sort(listEcc.begin(), listEcc.end());
  	vector<double>::iterator last = std::unique(listEcc.begin(), listEcc.end());
	listEcc.erase(last, listEcc.end());
}

double eccentricities(vector< vector <double> > &d, vector<double> &distribution, int i)
{
	/*
		eccentricities: 
		See "3.4.1 First Lower Bound" in "Quantitative comparison of protein isosurfaces 
 			with approximated Gromov-Wasserstein-distance" Definition 3.21

		Input 
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points
		distribution - vector containing the weight for the points
		i - integer in range 0 to number_of_selected_points, which is the
			index of the point for that the eccentricity is to calculate

		Output
		------
		Is a double that is the eccentricity of the i-th point

		Complexity
    		----------
		O(n), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > d{{0, 2, 2},
						{2, 0, 2},	
						{2, 2, 0}};
		vector <double> distribution{0.25,0.5,0.25};
		int i = 1;
		
		eccentricities(d, distribution, i);

		=> 1
	*/
	double s = 0;

	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		s += distribution[k] * d[i][k];
	}
	return(s);
}

void all_eccentricities(vector< vector <double> > &d, vector<double> &distribution,vector<double> &ecc)
{	
	/*
		all_eccentricities: 
		calculates all eccentricities for a given set of points

		Input 
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points
		distribution - vector containing the weight for the points
		ecc - is a vector, here all eccentrisities will be stored

		Output
		------
		No return

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > d{{0, 2, 2},
						{2, 0, 2},
						{2, 2, 0}};
		vector <double> distribution{0.25,0.5,0.25};
		vector <double> ecc(3);
		
		all_eccentricities(d, distribution, ecc);

		=> 1.5 1 1.5
	*/
	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		ecc[k] = eccentricities(d, distribution, k);
	}
}

void which_index(vector<double> &ecc, vector<double> &list, double u)
{
	/*
		which_index: 
		determines all indices of points that have an eccentricity lower than u

		Input 
		-----
		ecc - vector that contains all eccentricities of the points in an objects
		list - here the indices that have an eccentricity lower than u will be stored
		u - a double, that specifies how low the eccentricities must be

		Output
		------
		No return

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector <double> list;
		vector <double> ecc{1.5, 1, 1.5};
		double u = 1;

		which_index(ecc, list, u);
		
		=> 1
	*/
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
	/*
		distr_of_eccentricities: 
		Calculates the probability if one choose a point randomly, with the given distribution, 
		that the eccentricity is lower than u.

		Input 
		-----
		ecc - vector that contains all eccentricities of the points in an objects
		distribution - vector containing the weight for the points
		u - a double, that specifies how low the eccentricities must be

		Output
		------
		Is a double that is the distribution of eccentricities.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector <double> distribution{0.25,0.5,0.25};
		vector <double> ecc{1.5, 1, 1.5};
		double u = 1;
		
		distr_of_eccentricities(ecc, distribution, u);

		=> 0.5
	*/
	
	vector<double>list_index_of_points_ecc_le_u;
	double mu = 0;

	which_index(ecc, list_index_of_points_ecc_le_u, u);

	for(unsigned int i = 0; i < list_index_of_points_ecc_le_u.size(); i++)
	{
		mu = mu + distribution[list_index_of_points_ecc_le_u[i]];
	}
	return(mu);
}


double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, 
	   vector<double> &distribution2, int number_of_selected_points)
{
	/*
		flb: 
		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in 
			"Quantitative comparison of protein isosurfaces 
 			with approximated Gromov-Wasserstein-distance" from Felix Berens.

		Input 
		-----
		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		distribution1 - vector containing the weight for the points of the first object
		distribution2 - vector containing the weight for the points of the second object
		number_of_selected_points - a integer, for how many points the first lower bound will be calculated
		
		Output
		------
		Is a double that is the first lower bound.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > points1{{0,0,1},
						{1,0,0},
						{0,1,0}};
		vector < vector<double> > points2{{2,0,1},
						{1,0,2},
						{0,1,3}};
		vector <double> distribution1{0.4,0.2,0.2};
		vector <double> distribution2{0.25,0.5,0.25};
		int number_of_selected_points = 3;
		flb(points1, points2, distribution1, distribution2, number_of_selected_points);

		=> 0.372717
	*/
	vector < vector<double> > d1(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector <double> eccentricities1(number_of_selected_points);
	vector <double> eccentricities2(number_of_selected_points);
	vector<double> ecc1(number_of_selected_points);
	vector<double> ecc2(number_of_selected_points);
	vector<double> listEcc(number_of_selected_points+number_of_selected_points);
	double u;
	double out = 0;
	
	// calculation of the distance matrices
	dist_matrix(d1, points1);
	dist_matrix(d2, points2);

	// calculation of all eccentricities
	all_eccentricities(d1, distribution1,ecc1);
	all_eccentricities(d2, distribution2,ecc2);

	// mergin and sorting the list of eccentricities
	list_of_unique_eccent(ecc1, ecc2, listEcc);

	// calculating the lower bound
	for(unsigned int i = 0; i < listEcc.size() - 1; i++)
	{
		u = listEcc[i];
		out += (listEcc[i+1]-listEcc[i])*abs(distr_of_eccentricities(ecc1, distribution1, u) 
						     - distr_of_eccentricities(ecc2, distribution2, u));
	}
	return(0.5 * out);
}
#endif /* FIRST_LOWER_BOUND_H_ */
