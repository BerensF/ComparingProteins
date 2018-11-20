/*
* main.cpp
*
*  Created on: 07.05.2018
*      Author: Berens
*/

#include <iostream>
#include <vector>
#include <time.h>
#include <pthread.h>
#include "first_lower_bound.h"
#include "../reading.h"
#include "../random_points.h"

using namespace std;

#define NUM_THREADS 1

vector<double> out;
vector<double> distribution1;
vector<double> distribution2;
int rounds;
int number_of_selected_points;

struct thread_data {
	int  thread_id;
	vector< vector <double> > points1;
	vector< vector <double> > points2;
};
void *Thread(void *threadarg); // Threadfunction for multithreading

int main(int argc, char *argv[])
{
	string Dataname1 = argv[1];
	string Dataname2 = argv[2];
	number_of_selected_points = atoi(argv[3]);
	rounds = atoi(argv[4]);

	pthread_t threads[NUM_THREADS];
	struct thread_data td[NUM_THREADS];
	int rc;

	vector< vector <double> > points1(3);
	vector< vector <double> > points2(3);

	distribution1.resize(number_of_selected_points, 1.0/number_of_selected_points);
	distribution2.resize(number_of_selected_points, 1.0/number_of_selected_points);
	out.resize(rounds, -1);

	reading(Dataname1, points1); // reads the first file containing the 3D points
	if(points1[0].size() == 0){exit(-1);}
	reading(Dataname2, points2); // reads the second file containing the 3D points
	if(points2[0].size() == 0){exit(-1);}

	for(int i = 0; i < NUM_THREADS; i++ ) 
	{
		td[i].thread_id = i;
		td[i].points1 = points1;
		td[i].points2 = points2;

		rc = pthread_create(&threads[i], NULL, Thread, (void *)&td[i]);

		if (rc) {
			cerr << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_join(threads[i], NULL);
	}

	for(int i = 0; i< out.size();i++)
	{
		cout << out[i] << "\n";
	}
	return(0);
}


void *Thread(void *threadarg) {

	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;

	vector< vector <double> > selected_points1(3, vector<double>(number_of_selected_points));
	vector< vector <double> > selected_points2(3, vector<double>(number_of_selected_points));
	vector<unsigned int> random_indices_list;
	double SolutionOfBound;

	// /* initialize random seed: */
	srand(time(NULL) + my_data->thread_id);

	for(int j = 0; j < rounds/NUM_THREADS; j++)
	{		
		random_indices_list = random_points(my_data->points1[0].size(),number_of_selected_points); // Select number_of_selected_points indices randomly

		for(int i = 0; i < number_of_selected_points; i++)
		{
		selected_points1[0][i] = my_data->points1[0][random_indices_list[i]];
		selected_points1[1][i] = my_data->points1[1][random_indices_list[i]];
		selected_points1[2][i] = my_data->points1[2][random_indices_list[i]];
		}

		random_indices_list.clear();
		random_indices_list = random_points(my_data->points2[0].size(),number_of_selected_points); // Select number_of_selected_points indices randomly

		for(int i = 0; i < number_of_selected_points; i++)
		{
		selected_points2[0][i] = my_data->points2[0][random_indices_list[i]];
		selected_points2[1][i] = my_data->points2[1][random_indices_list[i]];
		selected_points2[2][i] = my_data->points2[2][random_indices_list[i]];
		}
		random_indices_list.clear();

		SolutionOfBound = flb(selected_points1, selected_points2, distribution1,distribution2, number_of_selected_points); // Calculates the first lower bound
		out[j + (my_data->thread_id)*(rounds)/NUM_THREADS] = SolutionOfBound;
	}
	pthread_exit(NULL);
	return(NULL);
}
