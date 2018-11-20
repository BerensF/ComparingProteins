/*
 * mathFunctions.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *
 *		Simple mathematics functions that are needed
 *
 */

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <cmath>

using namespace std;


double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2); // Calculates the euclidean distance for 2 3D points
void dist_matrix(vector < vector<double> > &d, vector< vector <double> > &points ); // Calculates a distance matrix for given 3D points


double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2)
{
	return(sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2)));
}

void dist_matrix(vector < vector<double> > &d, vector< vector <double> > &points )
{
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
