//
// Created by USER on 20/03/2024.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H


#include <vector>
#include <limits>
#include <iostream>
using namespace std;

class ZBuffer{
    double width;
    double height;

public:
    ZBuffer(int width, int height);

    std::vector<std::vector<double>> buffer;

    virtual ~ZBuffer();
};


#endif //ENGINE_ZBUFFER_H
