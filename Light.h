//
// Created by USER on 24/05/2024.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H
#include "Color.h"
#include "vector3d.h"
#include <iostream>
#include <list>

class Light {
public:
    Color ambientLight;
    Color diffuseLight;
    Color specularLight;
};



class PointLight: public Light{
public:
    Vector3D location;
    double spotAngle;
};

class InfLight: public Light{
public:
    Vector3D ldVector;
};

typedef std::list<Light> Lights3D;
#endif //ENGINE_LIGHT_H
