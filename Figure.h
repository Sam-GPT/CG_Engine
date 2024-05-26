//
// Created by USER on 6/03/2024.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H


#include "vector3d.h"
#include "Color.h"
#include "Face.h"


class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces ;
    Color ambientReflection;
    Color diffuseReflection;
    Color reflectionCoefficient;

};
typedef std::list<Figure>Figures3D;


#endif //ENGINE_FIGURE_H
