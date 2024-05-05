#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>

#include "cluster.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tbb.h"

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


        // This part should not be parallelized 
        // sum = tbb::parallel_reduce( 
        //     tbb::blocked_range<int>(0, total_values),
        //     0.0,
        //     [&]( const tbb::blocked_range<int> &r, double init) 
        //     -> double {
        //         for(int i = r.begin(); i != r.end(); ++i){
        //             init += pow(this->clusters[0].getCentralValue(i) - point.getValue(i), 2.0);
        //         }
        //         return init;
        //     },
        //     [](double x, double y) 
        //     -> double {
        //         return x + y;
        //     }
        // ); 

        for(int i = 0; i < total_values; i++)
        {
            sum += pow(clusters[0].getCentralValue(i) -
                       point.getValue(i), 2.0);
        }


        min_dist = sqrt(sum);

        for(int i = 1; i < K; i++)
        {
            double dist;

            // This part should not be parallelized 
            // sum = tbb::parallel_reduce( 
            //         tbb::blocked_range<int>(0, total_values),
            //         0.0,
            //         [&]( const tbb::blocked_range<int> &r, double init) 
            //         -> double {
            //             for(int j = r.begin(); j != r.end(); ++j){
            //                 init += pow(this->clusters[i].getCentralValue(j) - point.getValue(j), 2.0);
            //             }
            //             return init;
            //         },
            //         [](double x, double y) 
            //         -> double {
            //             return x + y;
            //         }
            //     );

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
        for(int i = 0; i < K; ++i)
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



            // // associates each point to the nearest center (sequential method 1)
            // for(int i = 0; i < total_points; i++)
            // {
            //     int id_old_cluster = points[i].getCluster();
            //     int id_nearest_center = getIDNearestCenter(points[i]);

            //     if(id_old_cluster != id_nearest_center)
            //     {
            //         if(id_old_cluster != -1)
            //             clusters[id_old_cluster].removePoint(points[i].getID());

            //         points[i].setCluster(id_nearest_center);
            //         clusters[id_nearest_center].addPoint(points[i]);
            //         done = false;
            //     }
            // }

        

            // // associates each point to the nearest center (sequential method 2)
            // std::vector<int> new_cluster_points(total_points, -2); 

            // for(int i = 0; i < total_points; ++i)
            // {
            //     int id_old_cluster = points[i].getCluster();
            //     int id_nearest_center = getIDNearestCenter(points[i]);

            //     if(id_old_cluster != id_nearest_center)
            //     {
            //         new_cluster_points.at(i) = id_nearest_center;
            //     }
            // }

            // for(int i = 0; i < total_points; ++i)
            // {
            //     if(new_cluster_points.at(i) != -2)
            //     {   
            //         int point_cluster = points[i].getCluster();

            //         if (point_cluster != -1)
            //             clusters[point_cluster].removePoint(points[i].getID());

            //         clusters[new_cluster_points.at(i)].addPoint(points[i]);
            //         points[i].setCluster(new_cluster_points.at(i));
            //         done = false;
            //     }

            // }



            // associates each point to the nearest center (parallel method 2)
            std::vector<int> new_cluster_points(total_points, -2); 

            tbb::parallel_for( tbb::blocked_range<int>(0, total_points),
                [&](const tbb::blocked_range<int> &r )
                {
                    for(int i = r.begin(); i < r.end(); ++i)
                    {
                        int id_old_cluster = points[i].getCluster();
                        int id_nearest_center = getIDNearestCenter(points[i]);

                        if(id_old_cluster != id_nearest_center)
                        {
                            new_cluster_points.at(i) = id_nearest_center;
                        }
                    }
                }

            );

            for(int i = 0; i < total_points; ++i)
            {
                if(new_cluster_points.at(i) != -2)
                {   
                    int point_cluster = points[i].getCluster();

                    if (point_cluster != -1)
                        clusters[point_cluster].removePoint(points[i].getID());

                    clusters[new_cluster_points.at(i)].addPoint(points[i]);
                    points[i].setCluster(new_cluster_points.at(i));
                    done = false;
                }

            }



            // for(int i = 0; i < total_points; ++i)
            // {
            //     Point& point_i = points[i];
            //     int point_cluster = point_i.getCluster();

            //     for(int c = 0; c < K; ++c)
            //     {
            //         if(new_cluster_points.at(i) != c)
            //         {   

            //             if (point_cluster != -1)
            //                 clusters[point_cluster].removePoint(point_i.getID());

            //             clusters[c].addPoint(point_i);
            //             point_i.setCluster(c);
            //             done = false;
            //         }
            //     }

            // }


            // // associates each point to the nearest center (parallel method 3)
            // for(int c = 0; c < K; ++c)
            // {
            //     for(int i = 0; i < total_points; ++i)
            //     {
            //         if(new_cluster_points.at(i) != c)
            //         {   
            //             int point_cluster = points[i].getCluster();

            //             if (point_cluster != -1)
            //                 clusters[point_cluster].removePoint(points[i].getID());

            //             clusters[c].addPoint(points[i]);
            //             points[i].setCluster(c);
            //             done = false;
            //         }

            //     }
            // }

            // for(int i = 0; i < total_points; ++i)
            // {
            //     Point point_i = points[i];
            //     if(new_cluster_points.at(i) != -2)
            //     {   
            //         int point_cluster = point_i.getCluster();

            //         if (point_cluster != -1)
            //             clusters[point_cluster].removePoint(point_i.getID());

            //         clusters[new_cluster_points.at(i)].addPoint(point_i);
            //         points[i].setCluster(new_cluster_points.at(i));
            //         done = false;
            //     }

            // }







            // // recalculating the center of each cluster (sequential method 1)
            // for(int i = 0; i < K; ++i)
            // {
            //     int total_points_cluster = clusters[i].getTotalPoints();

            //     if(total_points_cluster > 0)
            //     {   
            //         std::vector<double> sum(total_values, 0.0);
            //         for(int p = 0; p < total_points_cluster; ++p)
            //         {
            //             for(int j = 0; j < total_values; ++j)
            //             {
            //                 sum[j] += clusters[i].getPoint(p).getValue(j);
            //             }
            //         }

            //         for (int j = 0; j < total_values; ++j)
            //         {
            //             clusters[i].setCentralValue(j, sum[j] / total_points_cluster );
            //         }
                    
            //     }
                
            // }

            // recalculating the center of each cluster (parallel method 1)
            using vector_t = std::vector<double>;
            using private_vec_t = tbb::enumerable_thread_specific<vector_t>;


            for(int i = 0; i < K; i++)
            {
                int total_points_cluster = clusters[i].getTotalPoints();

                if(total_points_cluster > 0)
                { 
                    private_vec_t private_sum{total_values};

                    parallel_for(tbb::blocked_range<int>{0, total_points_cluster},
                        [&](const tbb::blocked_range<int> &r)
                        {
                            private_vec_t::reference my_sum = private_sum.local();

                            for(int p = r.begin(); p != r.end(); ++p)
                            {
                                for(int j = 0; j < total_values; ++j)
                                {
                                    my_sum[j] += this->clusters[i].getPoint(p).getValue(j);
                                }   
                            }
                        }
                        );

                    vector_t sum = private_sum.combine(
                        [] (vector_t a, vector_t b) -> vector_t 
                        {
                            std::transform
                            (
                                a.begin(), a.end(), b.begin(), a.begin(), std::plus<double> ()
                            );
                        return a;
                        }
                    );


                    for (int j = 0; j < total_values; ++j)
                    {
                        clusters[i].setCentralValue(j, sum[j] / total_points_cluster);
                    }

                    
                }
                
            }


            
            // // recalculating the center of each cluster (sequential method 2)
            // for(int i = 0; i < K; ++i)
            // {
            //     for(int j = 0; j < total_values; ++j)
            //     {
            //         int total_points_cluster = clusters[i].getTotalPoints();
            //         double sum = 0.0;

            //         if(total_points_cluster > 0)
            //         {
            //             for(int p = 0; p < total_points_cluster; ++p)
            //                 sum += clusters[i].getPoint(p).getValue(j);
            //             clusters[i].setCentralValue(j, sum / total_points_cluster);
            //         }
            //     }
            // }


            // // recalculating the center of each cluster  (parallel method 2)
            // for(int i = 0; i < K; ++i)
            // {
            //     for(int j = 0; j < total_values; ++j)
            //     {
            //         int total_points_cluster = this->clusters[i].getTotalPoints();

            //         if(total_points_cluster > 0)
            //         {
            //             double sum = tbb::parallel_reduce( 
            //                 tbb::blocked_range<int>(0, total_points_cluster),
            //                 0.0,
            //                 [&]( const tbb::blocked_range<int> &c, double init) -> double 
            //                 {
            //                     for(int p = c.begin(); p != c.end(); ++p)
            //                     {
            //                         init += this->clusters[i].getPoint(p).getValue(j);
            //                     }
            //                     return init;
            //                 },
            //                 [](double x, double y) -> double 
            //                 {
            //                     return x + y;
            //                 }
            //             ); 
                        
            //             this->clusters[i].setCentralValue(j, sum / total_points_cluster);
            //         }
            //     }
            // }



            // // recalculating the center of each cluster (parallel method 3)
            // using vector_t = std::vector<double>;
            // using private_vec_t = tbb::enumerable_thread_specific<vector_t>;
            
            // parallel_for(tbb::blocked_range<int>{0, K},
            //     [&](const tbb::blocked_range<int> &v)
            //     {
            //         for(int i = v.begin(); i != v.end(); ++i)
            //         {

            //             int total_points_cluster = clusters[i].getTotalPoints();

            //             if(total_points_cluster > 0)
            //             { 
            //                 private_vec_t private_sum{total_values};

            //                 parallel_for(tbb::blocked_range<int>{0, total_points_cluster},

            //                     [&](const tbb::blocked_range<int> &r)
            //                     {
            //                         private_vec_t::reference my_sum = private_sum.local();

            //                         for(int p = r.begin(); p != r.end(); ++p)
            //                         {
            //                             for(int j = 0; j < total_values; ++j)
            //                             {
            //                                 my_sum[j] += this->clusters[i].getPoint(p).getValue(j);
            //                             }   
            //                         }
            //                     }
            //                     );

            //                 vector_t sum = private_sum.combine(
            //                     [] (vector_t a, vector_t b) -> vector_t 
            //                     {
            //                         std::transform(a.begin(), a.end(), b.begin(), 
            //                             a.begin(), std::plus<double> ()
            //                             );
            //                     return a;
            //                     }
            //                 );


            //                 for (int j = 0; j < total_values; ++j)
            //                 {
            //                     clusters[i].setCentralValue(j, sum[j] / total_points_cluster);
            //                 }

            //             }
            //         }
                    
            //     }
            //     );


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
            //         cout << clusters[i].getPoint(j).getValue(p) << " ";

            //     string point_name = clusters[i].getPoint(j).getName();

            //     if(point_name != "")
            //         cout << "- " << point_name;

            //     cout << endl;
            // }

            // cout << "Cluster values: ";

            // for(int j = 0; j < total_values; j++)
            //     cout << clusters[i].getCentralValue(j) << " " ; 
            // cout << endl;           
            
        }

        cout << "TOTAL DITANCE OF ALL CLUSTERS = " << total_distance << endl;

        cout << "EXECUTION TIME PER ITERATION = "<<static_cast<float>(std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()/static_cast<float>(iter) )<<" (s)\n";


        cout << "TOTAL EXECUTION TIME = "<<std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()<<" (s)\n";
            
        cout << "TIME PHASE 1 = "<<std::chrono::duration_cast<std::chrono::seconds>(end_phase1-begin).count()<<" (s)\n";
        
        cout << "TIME PHASE 2 = "<<std::chrono::duration_cast<std::chrono::seconds>(end-end_phase1).count()<<" (s)\n";

    }
};