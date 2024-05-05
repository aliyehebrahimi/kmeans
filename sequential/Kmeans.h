#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>

#include "cluster.h"

using namespace std;

class KMeans
{
private:
    int K; // number of clusters
    int total_values, total_points, max_iterations;
    vector<Cluster> clusters;

    // return ID of nearest center (uses euclidean distance)
    int getIDNearestCenter(Point point)
    {
        double sum = 0.0, min_dist;
        int id_cluster_center = 0;

        for(int i = 0; i < total_values; i++)
        {
            sum += pow(clusters[0].getCentralValue(i) -
                       point.getValue(i), 2.0);
        }

        min_dist = sqrt(sum);

        for(int i = 1; i < K; i++)
        {
            double dist;
            sum = 0.0;

            for(int j = 0; j < total_values; j++)
            {
                sum += pow(clusters[i].getCentralValue(j) -
                           point.getValue(j), 2.0);
            }

            dist = sqrt(sum);

            if(dist < min_dist)
            {
                min_dist = dist;
                id_cluster_center = i;
            }
        }

        return id_cluster_center;
    }

public:
    KMeans(int K, int total_points, int total_values, int max_iterations)
    {
        this->K = K;
        this->total_points = total_points;
        this->total_values = total_values;
        this->max_iterations = max_iterations;
    }

    void run(vector<Point> & points)
    {
        auto begin = chrono::high_resolution_clock::now();
        
        if(K > total_points)
            return;

        vector<int> prohibited_indexes;

        // choose K distinct values for the centers of the clusters
        for(int i = 0; i < K; i++)
        {
            while(true)
            {
                int index_point = rand() % total_points;

                if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
                        index_point) == prohibited_indexes.end())
                {
                    prohibited_indexes.push_back(index_point);
                    points[index_point].setCluster(i);
                    Cluster cluster(i, points[index_point]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        auto end_phase1 = chrono::high_resolution_clock::now();
        
        int iter = 1;

        while(true)
        {
            bool done = true;

            // associates each point to the nearest center
            for(int i = 0; i < total_points; i++)
            {
                int id_old_cluster = points[i].getCluster();
                int id_nearest_center = getIDNearestCenter(points[i]);

                if(id_old_cluster != id_nearest_center)
                {
                    if(id_old_cluster != -1)
                        clusters[id_old_cluster].removePoint(points[i].getID());

                    points[i].setCluster(id_nearest_center);
                    clusters[id_nearest_center].addPoint(points[i]);
                    done = false;
                }
            }

            // recalculating the center of each cluster
            for(int i = 0; i < K; i++)
            {
                for(int j = 0; j < total_values; j++)
                {
                    int total_points_cluster = clusters[i].getTotalPoints();
                    double sum = 0.0;

                    if(total_points_cluster > 0)
                    {
                        for(int p = 0; p < total_points_cluster; p++)
                            sum += clusters[i].getPoint(p).getValue(j);
                        clusters[i].setCentralValue(j, sum / total_points_cluster);
                    }
                }
            }

            if(done == true || iter >= max_iterations)
            {
                cout << "Break in iteration " << iter << "\n\n";
                break;
            }

            iter++;
        }
        auto end = chrono::high_resolution_clock::now();

        // shows elements of clusters
        double total_distance = 0.0;

        for(int i = 0; i < K; i++)
        {
            int total_points_cluster =  clusters[i].getTotalPoints();
            double distance = clusters[i].totalDistance();
            total_distance += distance;

            // cout << "Cluster " << clusters[i].getID() + 1 << "\n";
            // cout << "TOTAL DITANCE = " << distance << "\n\n";

            // for(int j = 0; j < total_points_cluster; j++)
            // {
            //     cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
            //     for(int p = 0; p < total_values; p++)
            //         cout << "(" << p<< ")"<<clusters[i].getPoint(j).getValue(p) << ", ";

            //     string point_name = clusters[i].getPoint(j).getName();

            //     if(point_name != "")
            //         cout << "- " << point_name;

            //     cout << endl;
            // }

            // cout << "Cluster values: ";

            // for(int j = 0; j < total_values; j++)
            //     cout << clusters[i].getCentralValue(j) << ", " ; 
            // cout << endl;           
            
        }

        cout << "TOTAL DITANCE OF ALL CLUSTERS = " << total_distance << endl;

        cout << "EXECUTION TIME PER ITERATION = "<<static_cast<float>(std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()/static_cast<float>(iter) )<<" (s)\n";


        cout << "TOTAL EXECUTION TIME = "<<std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()<<" (s)\n";
            
        cout << "TIME PHASE 1 = "<<std::chrono::duration_cast<std::chrono::seconds>(end_phase1-begin).count()<<" (s)\n";
        
        cout << "TIME PHASE 2 = "<<std::chrono::duration_cast<std::chrono::seconds>(end-end_phase1).count()<<" (s)\n";

    }
};