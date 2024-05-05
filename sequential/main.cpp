// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>

#include "Kmeans.h"

using namespace std;


int main(int argc, char *argv[])
{
	srand (time(NULL));

	int total_points, total_values, K, max_iterations, has_name, thread_num;

	cin >> total_points >> total_values >> K >> max_iterations >> has_name >> thread_num;

	vector<Point> points;
	string point_name;

	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;

		for(int j = 0; j < total_values; j++)
		{
			double value;
			cin >> value;
			values.push_back(static_cast<double>(value) );
		}

		if(has_name)
		{
			cin >> point_name;
			Point p(i, values, point_name);
			points.push_back(p);
		}
		else
		{
			Point p(i, values);
			points.push_back(p);
		}
	}

	cout << "TOTAL NUMBER OF CLUSTERS: " << K << endl;
	cout << "TOTAL NUMBER OF POINTS: " << total_points << endl;
    cout << "TOTAL NUMBER OF FEATURES: " << total_values << endl << endl;

    cout << "MAX NUMBER OF ITERATIONS: " << max_iterations << endl << endl;

	KMeans kmeans(K, total_points, total_values, max_iterations);
	kmeans.run(points);

	return 0;
}
