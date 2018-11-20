/*
 * Reading.h
 *
 *  Created on: 25.03.2018
 *      Author: Berens
 *
 *		This function reads all points of a given datafile [Dataname] and stores them to [points].
 *		The datafile has to start with a one line header. In the following lines are the points, one line for one point.
 *		The coordinates of the points have to be seperated by points.
 *		Example:
 *		"x";"y";"z"
 *		66;71.9759158681157;40
 *		66;72;39.9936496036735
 *		65.8866023427208;72;40
 *		66.0721155206865;72;40
 *		65;72.3806364831402;40
 *		65;73;39.920087406864
 *		....
 */

#ifndef READING_H_
#define READING_H_

#include "string.h"
#include "stdlib.h"
#include "fstream"

using namespace std;

void reading(string Dataname, vector< vector <double> > &points)
{
	ifstream input(Dataname.c_str());
	string word;
	string line;
	size_t pos_sep; // Position of separator
//	int number_of_lines;
	char zwsp[256];

	if (!input)
	{
		cerr << "Fehler beim Oeffnen der Datei " << Dataname << "\n";
	}else{
		input >> word; // header

//		cout << "Starting reading points from " << Dataname << "\n";

		while(!input.eof())
		{
			input >> line;

//			x value
			pos_sep = line.find(";");
			word = line.substr (0,pos_sep);
			strcpy(zwsp, word.c_str());
			points[0].push_back(atof(zwsp));
			line.erase(0,pos_sep+1);

//			 y value
			pos_sep = line.find(";");
			word = line.substr (0,pos_sep);
			strcpy(zwsp, word.c_str());
			points[1].push_back(atof(zwsp));
			line.erase(0,pos_sep+1);

//			z value
			pos_sep = line.find(";");
			word = line.substr (0,pos_sep);
			strcpy(zwsp, word.c_str());
			points[2].push_back(atof(zwsp));
			line.erase(0,pos_sep+1);
		}

//		cout << "End of reading points.\n";
	}
}


#endif /* READING_H_ */
