/*
 * mathFunctions.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *
 *		Simple mathematics functions that are needed
 */

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <cmath>

using namespace std;


double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2);
void dist_matrix(vector < vector<double> > &d, vector< vector <double> > &points );


double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2)
{
	/*
		euclidean_distance:
		Calculates the euclidean distance for 2 3D points

		Input 
		-----
		x1, y1, z1 - are the coordinates of the first point, all types are double
		x2, y2, z2 - are the coordinates of the second point, all types are double

		Output
		------
		The output is a double, which is the euclidean distance between (x1,y1,z1) and (x2,y2,z2)
		
		Complexity
    		----------
		O(1)

		Example
		-------
		double x1 = 0;
		double y1 = 0;
		double z1 = 0;
		double x2 = 1;
		double y2 = 1;
		double z2 = 1;
		euclidean_distance(x1,y1,z1,x2,y2,z2);

		=> 1.73205 
	*/
	return(sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2)));
}

void dist_matrix(vector < vector<double> > &d, vector< vector <double> > &points )
{
	/*
		dist_matrix:
		Calculates a distance matrix for given 3D points

		Input 
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points.
			The distance matrix will be stored there.
		points - is from type vector, with dimenson number_of_selected_points x 3.
			The points for which the distance matrix will be calculated are stored here.

		Output
		------
		No return
		
		Complexity
    		----------
		O(number_of_selected_points)

		Example
		-------
		vector < vector<double> > d(3, vector<double>(3,0));
		vector < vector<double> > points{ { 1, 0, 0 }, 
						  { 0, 0, 1 }, 		
						  { 0, 1, 0 } }; 
		dist_matrix(d,points);

			(0	 1.41421 1.41421)
		=> d =	(1.41421 0	 1.41421)
			(1.41421 1.41421 0	)
	*/
	for(unsigned int i=0; i < d.size(); i++)
	{
		for(unsigned int j=i;j< d.size(); j++)
		{
			d[i][j] = euclidean_distance(points[0][i],points[1][i],points[2][i],points[0][j],points[1][j],points[2][j]);
			d[j][i] = d[i][j];
		}
	}
}


#endif /* MATHFUNCTIONS_H_ */
