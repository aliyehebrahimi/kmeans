#include <iostream>
#include <vector>
#include <stdlib.h>

#include "point.h"

using namespace std;

class Cluster
{
private:
    int id_cluster;
    vector<double> central_values;
    vector<Point> points;

public:
    Cluster(int id_cluster, Point point)
    {
        this->id_cluster = id_cluster;

        int total_values = point.getTotalValues();

        for(int i = 0; i < total_values; i++)
            central_values.push_back(point.getValue(i));

        points.push_back(point);
    }

    void addPoint(Point point)
    {
        points.push_back(point);
    }

    bool removePoint(int id_point)
    {
        int total_points = points.size();

        for(int i = 0; i < total_points; i++)
        {
            if(points[i].getID() == id_point)
            {
                points.erase(points.begin() + i);
                return true;
            }
        }
        return false;
    }

    double getCentralValue(int index)
    {
        return central_values[index];
    }

    void setCentralValue(int index, double value)
    {
        central_values[index] = value;
    }

    Point getPoint(int index)
    {
        return points[index];
    }

    int getTotalPoints()
    {
        return points.size();
    }

    int getID()
    {
        return id_cluster;
    }

    double totalDistance()
    {
            if (this->points.empty())
            {
                return 0.0;
            }

            double sum = 0.0;
            int total_points = this->points.size();
            int total_values = this->points.at(0).getTotalValues();


            double total_sum = 0;

            for(int i = 0; i < total_points; ++i)
            {
                double point_sum = 0;
                for(int j = 0; j < total_values; ++j)
                {
                    point_sum += pow(this->central_values[j] - 
                                    this->points.at(i).getValue(j), 2.0);
                }
                total_sum += point_sum;
            }

        return total_sum;
    }

};