#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include "l_parser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <cmath>
#include <set>
#include <stack>
#include "vector3d.h"
#include "Figure.h"
#include <string>
#include <sstream>
#include "ZBuffer.h"
#include <algorithm>
#define PI 3.14159265358979323846


using Lines2D = std::list<Line2D>;
void getXmin(const Line2D& line, double& cur_min){
    if(line.p1.x < cur_min){
        cur_min = line.p1.x;
    }
    else if(line.p2.x < cur_min){
        cur_min = line.p2.x;
    }
}

void getXmax(const Line2D& line, double& cur_max){
    if(line.p1.x > cur_max){
        cur_max = line.p1.x;
    }
    else if(line.p2.x > cur_max){
        cur_max= line.p2.x;
    }
}

void getYmin(const Line2D& line, double& cur_min){
    if(line.p1.y < cur_min){
        cur_min = line.p1.y;
    }
    else if(line.p2.y < cur_min){
        cur_min = line.p2.y;
    }
}

void getYmax(const Line2D& line, double& cur_max){
    if(line.p1.y > cur_max){
        cur_max = line.p1.y;
    }
    else if(line.p2.y > cur_max){
        cur_max= line.p2.y;
    }
}
void mult_cords_with_SF(double& d, Lines2D& lines){
    for(auto& line : lines){
        line.p1.x *=d;
        line.p1.y *= d;
        line.p2.x *=d;
        line.p2.y *= d;
    }
}

std::string Replace(std::string& string, const LParser::LSystem2D &l_system){
    std::string newstring = "";

    for(auto& i : string){
        if(i== '+' || i == '-' || i == '(' || i == ')' || i == '^' || i == '&' || i == '\'' || i == '/'){
            newstring+=i;
        }
        else{
            newstring+= l_system.get_replacement(i);
        }

    }

    return newstring;

}
std::string FullString(const LParser::LSystem2D &l_system){
    std::string cur_string = l_system.get_initiator();

    for(size_t i = 0; i < l_system.get_nr_iterations(); i++){
        std::string newstring = Replace(cur_string, l_system);
        cur_string = newstring;
    }

    return cur_string;
}

std::string Replace(std::string& string, const LParser::LSystem3D &l_system){
    std::string newstring;
    for(auto& i : string){
        if(i== '+' || i == '-' || i == '(' || i == ')' || i == '^' || i == '&' || i == '\\' || i == '/' || i == '|'){
            newstring+=i;
        }
        else{
            newstring+= l_system.get_replacement(i);
            
        }

    }
    return newstring;

}
std::string FullString(const LParser::LSystem3D &l_system){
    std::string cur_string = l_system.get_initiator();

    for(size_t i = 0; i < l_system.get_nr_iterations(); i++){
        std::string newstring = Replace(cur_string, l_system);
        cur_string = newstring;

    }
    return cur_string;
}

double degreesToRadians(const double& degrees) {
    return degrees * (PI / 180.0);
}

void draw_zbuf_line(ZBuffer &buffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const img::Color& color){
    if (x0 >= image.get_width() || y0 >= image.get_height() || x1 >= image.get_width() || y1 > image.get_height()) {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << image.get_width() << " and height " << image.get_height();
        throw std::runtime_error(ss.str());
    }

    if (x0 == x1)
    {
        //special case for x0 == x1
        unsigned int aantal_pixels = std::max(y0 ,y1) - std::min(y0, y1) + 1;
        unsigned int current_pixelNum = aantal_pixels;

        if(y0 > y1){
            std::swap(z0, z1);
        }

        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++){
            double p = (double)current_pixelNum/aantal_pixels;
            double zb_waarde = (p/z0) + ((1-p)/z1);


            if(zb_waarde < buffer.buffer[x0][i]){
                (image)(x0, i) = color;
                buffer.buffer[x0][i] = zb_waarde;

            }
            current_pixelNum--;

        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        unsigned int aantal_pixels = std::max(x0 ,x1) -std::min(x0, x1) + 1;
        unsigned int current_pixelNum = aantal_pixels;

        if(x0 > x1){
            std::swap(z0, z1);
        }

        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++){
            double p = (double)current_pixelNum/aantal_pixels;
            double zb_waarde = (p/z0) + ((1-p)/z1);
            if(zb_waarde < buffer.buffer[i][y0]){
                (image)(i, y0) = color;
                buffer.buffer[i][y0] = zb_waarde;

            }
            current_pixelNum--;
        }


    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            unsigned int aantal_pixels = std::max(x1, x0) - std::min(x1, x0) + 1;
            unsigned int current_pixelNum = aantal_pixels;
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double p = (double)current_pixelNum/aantal_pixels;
                double zb_waarde = (p/z0) + ((1-p)/z1);
                if(zb_waarde < buffer.buffer[x0 + i][(unsigned int) round(y0 + m * i)]){
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                    buffer.buffer[x0 + i][(unsigned int) round(y0 + m * i)] = zb_waarde;

                }
                current_pixelNum--;
            }

        }
        else if (m > 1.0)
        {
            unsigned int aantal_pixels = std::max(y1, y0) - std::min(y1, y0) + 1;
            unsigned int current_pixelNum = aantal_pixels;


            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double p = (double)current_pixelNum/aantal_pixels;
                double zb_waarde = (p/z0) + ((1-p)/z1);
                if(zb_waarde < buffer.buffer[(unsigned int) round(x0 + (i / m))][y0 + i]){
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                    buffer.buffer[(unsigned int) round(x0 + (i / m))][y0 + i] = zb_waarde;

                }
                current_pixelNum--;
            }

        }
        else if (m < -1.0)
        {
            unsigned int aantal_pixels = std::max(y0, y1) - std::min(y0, y1) + 1;
            unsigned int current_pixelNum = aantal_pixels;


            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double p = (double)current_pixelNum/aantal_pixels;
                double zb_waarde = (p/z0) + ((1-p)/z1);
                if(zb_waarde < buffer.buffer[(unsigned int) round(x0 - (i / m))][y0 - i]){
                    (image)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                    buffer.buffer[(unsigned int) round(x0 - (i / m))][y0 - i] = zb_waarde;
                }
                current_pixelNum--;
            }

        }
    }
}

void Calculate(double& imageX, double& imageY, double& dx, double& dy, double& schaalfactor, const Lines2D &lines,const int size, const ini::Configuration &configuration){
    double xmin = lines.begin()->p1.x;
    double xmax = lines.begin()->p1.x;
    double ymin = lines.begin()->p1.y;
    double ymax = lines.begin()->p1.y;
    for(auto& line : lines){
        getXmin(line, xmin);
        getXmax(line, xmax);
        getYmin(line, ymin);
        getYmax(line, ymax);
    }
    double xrange = xmax - xmin;
    double yrange = ymax - ymin;


    imageX = size * (xrange/ fmax(xrange,yrange));
    imageY = size * (yrange/ fmax(xrange,yrange));


    schaalfactor = 0.95 * (imageX/xrange);


    double DCx = schaalfactor * ((xmin + xmax)/2);
    double DCy = schaalfactor * ((ymin + ymax)/2);

    dx = (imageX/2) - DCx;
    dy = (imageY/2) - DCy;


}

vector<Face> triangulate(const Face& face){
    vector<Face> faces;

    for(size_t i = 1 ; i <= face.point_indexes.size()-2 ; i++){
        Face newFace;
        newFace.point_indexes.push_back(face.point_indexes[0]);
        newFace.point_indexes.push_back(face.point_indexes[i]);
        newFace.point_indexes.push_back(face.point_indexes[i + 1]);

        faces.push_back(newFace);

    }
    return faces;
}

img::EasyImage draw2DLines(const Lines2D &lines,const int size, const ini::Configuration &configuration){
    double imageX;
    double imageY;
    double dx;
    double dy;
    double schaalfactor;
    Lines2D newLines;

    double xmin = lines.begin()->p1.x;
    double xmax = lines.begin()->p1.x;
    double ymin = lines.begin()->p1.y;
    double ymax = lines.begin()->p1.y;
    for(auto& line : lines){
        getXmin(line, xmin);
        getXmax(line, xmax);
        getYmin(line, ymin);
        getYmax(line, ymax);
    }

    double xrange = xmax - xmin;
    double yrange = ymax - ymin;


    imageX = size * (xrange/ fmax(xrange,yrange));
    imageY = size * (yrange/ fmax(xrange,yrange));


    schaalfactor = 0.95 * (imageX/xrange);


    double DCx = schaalfactor * ((xmin + xmax)/2);
    double DCy = schaalfactor * ((ymin + ymax)/2);

    dx = (imageX/2) - DCx;
    dy = (imageY/2) - DCy;



    newLines = lines;

    mult_cords_with_SF(schaalfactor, newLines);


    for(auto &line : newLines){
        line.p1.x+=dx;
        line.p1.x= ceil(line.p1.x);

        line.p2.x+=dx;
        line.p2.x=ceil(line.p2.x);


        line.p1.y+=dy;
        line.p1.y=ceil(line.p1.y);

        line.p2.y+=dy;
        line.p2.y= ceil(line.p2.y);

    }

    img::EasyImage image(imageX, imageY);
    std::vector<double> bgcolor = configuration["General"]["backgroundcolor"];
    image.clear(img::Color(bgcolor[0]*255, bgcolor[1]*255, bgcolor[2]*255));
    ZBuffer buffer(imageX, imageY);



    for(auto &line : newLines){
        if(configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){
            draw_zbuf_line(buffer, image, line.p1.x, line.p1.y, line.z1, line.p2.x, line.p2.y , line.z2, img::Color((line.color.red)*255,(line.color.green)*255, (line.color.blue)*255));
        }
        else{
            image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y,img::Color((line.color.red)*255,(line.color.green)*255, (line.color.blue)*255));
        }


    }
    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system , const ini::Configuration &configuration){
    Lines2D lines;
    std::stack<std::pair<Point2D, double>> myStack;

    double cur_angle = l_system.get_starting_angle();

    double angle = l_system.get_angle();

    std::string strings = FullString(l_system);

    std::set<char> alfabet = l_system.get_alphabet();

    Point2D cur_position;
    cur_position.x = 0;
    cur_position.y = 0;


    for(auto& i : strings){
        if(alfabet.find(i)!= alfabet.end() && l_system.draw(i)){
            Point2D new_position;
            new_position.x = cur_position.x + cos(degreesToRadians(cur_angle));
            new_position.y = cur_position.y + sin(degreesToRadians(cur_angle));

            Line2D newLine;
            newLine.p1 = cur_position;
            newLine.p2 = new_position;

            std::vector<double>vec_color =  configuration["2DLSystem"]["color"].as_double_tuple_or_die();

            newLine.color.red = vec_color[0];
            newLine.color.green = vec_color[1];
            newLine.color.blue = vec_color[2];

            cur_position.x = new_position.x;
            cur_position.y = new_position.y;
            lines.push_back(newLine);



        }
        else if(i == '+'){
            cur_angle+= angle;
        }
        else if(i == '-'){
            cur_angle-=angle;
        }
        else if(i == '('){
            std::pair<Point2D, double> topush = {cur_position, cur_angle};
            myStack.push(topush);
        }
        else{
            std::pair<Point2D, double> poped = myStack.top();
            myStack.pop();

            cur_position.x = poped.first.x;
            cur_position.y = poped.first.y;

            cur_angle = poped.second;
        }
    }


    return lines;

}

Figure draw3DLSystem(const LParser::LSystem3D &l_system , const ini::Configuration &configuration, std::string& fignum){
    Figure figuur;

    figuur.color.red = configuration[fignum]["color"].as_double_tuple_or_die()[0];
    figuur.color.green = configuration[fignum]["color"].as_double_tuple_or_die()[1];
    figuur.color.blue = configuration[fignum]["color"].as_double_tuple_or_die()[2];

    std::stack<std::vector<Vector3D>> myStack;

    double delta = l_system.get_angle();

    double delta_radians = degreesToRadians(delta);

    std::string strings = FullString(l_system);

    std::set<char> alfabet = l_system.get_alphabet();

    Vector3D cur_position = Vector3D::point(0,0,0);

    Vector3D vec_H = Vector3D::vector(1,0,0);

    Vector3D vec_L = Vector3D::vector(0,1,0);

    Vector3D vec_U = Vector3D::vector(0,0,1);


    for(auto& i : strings){
        Vector3D oud_hvector = vec_H;
        Vector3D oud_lvector = vec_L;
        Vector3D oud_uvector = vec_U;
        if(alfabet.find(i)!= alfabet.end() && l_system.draw(i)){

            figuur.points.push_back(cur_position);

            Vector3D new_point = cur_position + vec_H;

            figuur.points.push_back(new_point);

            cur_position = new_point;

            Face new_face;
            new_face.point_indexes.push_back(figuur.points.size()-2);
            new_face.point_indexes.push_back(figuur.points.size()-1);
            figuur.faces.push_back(new_face);


        }
        else if(i == '+'){
            vec_H = oud_hvector*cos(delta_radians) + oud_lvector*sin(delta_radians);
            vec_L = -oud_hvector*sin(delta_radians) + oud_lvector*cos(delta_radians);
        }
        else if(i =='-'){
            vec_H = oud_hvector*cos(-delta_radians) + oud_lvector*sin(-delta_radians);
            vec_L = -oud_hvector*sin(-delta_radians) + oud_lvector*cos(-delta_radians);
        }
        else if(i == '^'){
            vec_H = oud_hvector*cos(delta_radians) + oud_uvector*sin(delta_radians);
            vec_U = -oud_hvector*sin(delta_radians) + oud_uvector*cos(delta_radians);
        }
        else if(i == '&'){
            vec_H = oud_hvector*cos(-delta_radians) + oud_uvector*sin(-delta_radians);
            vec_U = -oud_hvector*sin(-delta_radians) + oud_uvector*cos(-delta_radians);
        }
        else if(i == '\\'){
            vec_L = oud_lvector*cos(delta_radians) - oud_uvector* sin(delta_radians);
            vec_U = oud_lvector*sin(delta_radians) + oud_uvector* cos(delta_radians);
        }
        else if(i == '/'){
            vec_L = oud_lvector*cos(-delta_radians) - oud_uvector* sin(-delta_radians);
            vec_U = oud_lvector*sin(-delta_radians) + oud_uvector* cos(-delta_radians);
        }
        else if(i == '|'){
            vec_H *= -1;
            vec_L *= -1;
        }
        else if(i == '('){
            std::vector<Vector3D> topush = {cur_position, vec_H, vec_L, vec_U};
            myStack.push(topush);
        }
        else{
            std::vector<Vector3D> poped = myStack.top();
            myStack.pop();

            cur_position = poped[0];
            vec_H = poped[1];
            vec_L = poped[2];
            vec_U = poped[3];

        }
    }

    return figuur;
}




Matrix scaleFigure(const double scale){
    Matrix m;

    m(1,1) = scale;
    m(2,2) = scale;
    m(3,3) = scale;
    return m;
}

Matrix rotateX(const double angle){
    Matrix m;

    m(2,2) = cos(angle);
    m(2,3) = sin(angle);
    m(3,2) = -1 * sin(angle);
    m(3,3) = cos(angle);

    return m;
}

Matrix rotateY(const double angle){
    Matrix m;

    m(1,1) = cos(angle);
    m(1, 3) = -1 * sin(angle);
    m(3,1) = sin(angle);
    m(3,3) = cos(angle);

    return m;
}

Matrix rotateZ(const double angle){
    Matrix m;

    m(1,1) = cos(angle);
    m(1,2) = sin(angle);
    m(2,1) = -1 * sin(angle);
    m(2,2) = cos(angle);

    return m;
}

Matrix translate(const Vector3D &vector){
    Matrix m;
    m(4,1) = vector.x;
    m(4,2) = vector.y;
    m(4,3) = vector.z;

    return m;
}

void applyTransformation(Figure &fig, const Matrix &m){
    for(auto &point : fig.points){
        point*=m;
    }

}

Point2D doProjection(const Vector3D &point,const double d){
    Point2D punt{};

    punt.x = (d*point.x)/-point.z;
    punt.y = (d*point.y)/-point.z;

    return punt;
}

Lines2D doProjection(const Figures3D &figs) {
    Lines2D lines;
    int fig_num = 0;
    for (const auto &fig : figs) {
        std::string fignum = "Figure"+ std::to_string(fig_num);
        for (const auto &vlak : fig.faces) {
            for (size_t i = 0; i < vlak.point_indexes.size() - 1; ++i) {
                int begin = vlak.point_indexes[i];
                int end = vlak.point_indexes[i + 1];


                Line2D newline{};
                newline.p1 = doProjection(fig.points[begin], 1);
                newline.p2 = doProjection(fig.points[end], 1);

                newline.z1 = fig.points[begin].z;
                newline.z2 = fig.points[end].z;

                newline.color.red = fig.color.red;
                newline.color.green = fig.color.green;
                newline.color.blue = fig.color.blue;


                lines.push_back(newline);
            }


            int begin = vlak.point_indexes.back();
            int end = vlak.point_indexes.front();


            Line2D newline{};
            newline.p1 = doProjection(fig.points[begin], 1);
            newline.p2 = doProjection(fig.points[end], 1);
            newline.z1 = fig.points[begin].z;
            newline.z2 = fig.points[end].z;

            newline.color.red = fig.color.red;
            newline.color.green = fig.color.green;
            newline.color.blue = fig.color.blue;


            lines.push_back(newline);

        }
        fig_num++;

    }

    return lines;
}

void toPolar(const Vector3D &point,double &theta,double &phi,double &r){
    r = sqrt(pow(point.x, 2)+ pow(point.y, 2)+ pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos(point.z/r);
}

Matrix eyePointTrans(const Vector3D &eyepoint){
    double rpoint;
    double theta;
    double phi;

    toPolar(eyepoint, theta, phi, rpoint);

    Matrix m;
    m(1,1) = -1 * sin(theta);
    m(1,2) = -1* cos(theta)* cos(phi);
    m(1,3) = cos(theta)*sin(phi);

    m(2,1) = cos(theta);
    m(2,2) = -1* sin(theta)* cos(phi);
    m(2,3) = sin(theta)*sin(phi);

    m(3,2) = sin(phi);
    m(3,3) = cos(phi);

    m(4,3) = -rpoint;

    return m;

}



Figure createTetrahedron(Color kleur){
    Figure t_Figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Vector3D p1 = Vector3D::point(1, -1, -1);
    Vector3D p2 = Vector3D::point(-1, 1, -1);
    Vector3D p3 = Vector3D::point(1, 1, 1);
    Vector3D p4 = Vector3D::point(-1, -1, 1);

    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);

    Face face1;
    face1.point_indexes.push_back(0);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(1);

    Face face2;

    face2.point_indexes.push_back(1);
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(3);

    Face face3;

    face3.point_indexes.push_back(0);
    face3.point_indexes.push_back(1);
    face3.point_indexes.push_back(3);

    Face face4;

    face4.point_indexes.push_back(0);
    face4.point_indexes.push_back(2);
    face4.point_indexes.push_back(3);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);
    faces.push_back(face4);

    t_Figure.points = points;
    t_Figure.faces = faces;
    t_Figure.color = kleur;

    return t_Figure;
}

Figure createIcosahedron(Color kleur){
    Figure I_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Vector3D p1 = Vector3D::point(0, 0, sqrt(5)/2);
    points.push_back(p1);
    for(int i = 2; i < 7 ;i++){
        Vector3D p = Vector3D::point(cos(((i - 2)*(2*PI))/5), sin(((i - 2)*(2*PI))/5), 0.5);
        points.push_back(p);
    }

    for(int i = 7; i < 12; i++){
        Vector3D p = Vector3D::point(cos((PI/5) + ((i-7)*(2*PI))/5), sin((PI/5) + ((i-7)*(2*PI))/5), -0.5);
        points.push_back(p);
    }

    Vector3D p12 = Vector3D::point(0, 0, -sqrt(5)/2);
    points.push_back(p12);

    Face face1;
    face1.point_indexes.push_back(0);
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);

    Face face2;
    face2.point_indexes.push_back(0);
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(3);

    Face face3;
    face3.point_indexes.push_back(0);
    face3.point_indexes.push_back(3);
    face3.point_indexes.push_back(4);

    Face face4;
    face4.point_indexes.push_back(0);
    face4.point_indexes.push_back(4);
    face4.point_indexes.push_back(5);

    Face face5;
    face5.point_indexes.push_back(0);
    face5.point_indexes.push_back(5);
    face5.point_indexes.push_back(1);

    Face face6;
    face6.point_indexes.push_back(1);
    face6.point_indexes.push_back(6);
    face6.point_indexes.push_back(2);

    Face face7;
    face7.point_indexes.push_back(2);
    face7.point_indexes.push_back(6);
    face7.point_indexes.push_back(7);

    Face face8;
    face8.point_indexes.push_back(2);
    face8.point_indexes.push_back(7);
    face8.point_indexes.push_back(3);

    Face face9;
    face9.point_indexes.push_back(3);
    face9.point_indexes.push_back(7);
    face9.point_indexes.push_back(8);

    Face face10;
    face10.point_indexes.push_back(3);
    face10.point_indexes.push_back(8);
    face10.point_indexes.push_back(4);

    Face face11;
    face11.point_indexes.push_back(4);
    face11.point_indexes.push_back(8);
    face11.point_indexes.push_back(9);

    Face face12;
    face12.point_indexes.push_back(4);
    face12.point_indexes.push_back(9);
    face12.point_indexes.push_back(5);

    Face face13;
    face13.point_indexes.push_back(5);
    face13.point_indexes.push_back(9);
    face13.point_indexes.push_back(10);

    Face face14;
    face14.point_indexes.push_back(5);
    face14.point_indexes.push_back(10);
    face14.point_indexes.push_back(1);

    Face face15;
    face15.point_indexes.push_back(1);
    face15.point_indexes.push_back(10);
    face15.point_indexes.push_back(6);

    Face face16;
    face16.point_indexes.push_back(11);
    face16.point_indexes.push_back(7);
    face16.point_indexes.push_back(6);

    Face face17;
    face17.point_indexes.push_back(11);
    face17.point_indexes.push_back(8);
    face17.point_indexes.push_back(7);

    Face face18;
    face18.point_indexes.push_back(11);
    face18.point_indexes.push_back(9);
    face18.point_indexes.push_back(8);

    Face face19;
    face19.point_indexes.push_back(11);
    face19.point_indexes.push_back(10);
    face19.point_indexes.push_back(9);

    Face face20;
    face20.point_indexes.push_back(11);
    face20.point_indexes.push_back(6);
    face20.point_indexes.push_back(10);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);
    faces.push_back(face4);
    faces.push_back(face5);
    faces.push_back(face6);
    faces.push_back(face7);
    faces.push_back(face8);
    faces.push_back(face9);
    faces.push_back(face10);
    faces.push_back(face11);
    faces.push_back(face12);
    faces.push_back(face13);
    faces.push_back(face14);
    faces.push_back(face15);
    faces.push_back(face16);
    faces.push_back(face17);
    faces.push_back(face18);
    faces.push_back(face19);
    faces.push_back(face20);

    I_figure.points = points;
    I_figure.faces = faces;
    I_figure.color = kleur;

    return I_figure;

}

Figure creatOctahedron(Color kleur){
    Figure o_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Vector3D p1 = Vector3D::point(1, 0, 0);
    Vector3D p2 = Vector3D::point(0, 1, 0);
    Vector3D p3 = Vector3D::point(-1, 0, 0);
    Vector3D p4 = Vector3D::point(0, -1, 0);
    Vector3D p5 = Vector3D::point(0, 0, -1);
    Vector3D p6 = Vector3D::point(0, 0, 1);

    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);
    points.push_back(p5);
    points.push_back(p6);

    Face face1;
    face1.point_indexes.push_back(0);
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(5);

    Face face2;
    face2.point_indexes.push_back(1);
    face2.point_indexes.push_back(2);
    face2.point_indexes.push_back(5);

    Face face3;
    face3.point_indexes.push_back(2);
    face3.point_indexes.push_back(3);
    face3.point_indexes.push_back(5);

    Face face4;
    face4.point_indexes.push_back(3);
    face4.point_indexes.push_back(0);
    face4.point_indexes.push_back(5);

    Face face5;
    face5.point_indexes.push_back(1);
    face5.point_indexes.push_back(0);
    face5.point_indexes.push_back(4);

    Face face6;
    face6.point_indexes.push_back(2);
    face6.point_indexes.push_back(1);
    face6.point_indexes.push_back(4);

    Face face7;
    face7.point_indexes.push_back(3);
    face7.point_indexes.push_back(2);
    face7.point_indexes.push_back(4);

    Face face8;
    face8.point_indexes.push_back(0);
    face8.point_indexes.push_back(3);
    face8.point_indexes.push_back(4);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);
    faces.push_back(face4);
    faces.push_back(face5);
    faces.push_back(face6);
    faces.push_back(face7);
    faces.push_back(face8);

    o_figure.points = points;
    o_figure.faces = faces;
    o_figure.color = kleur;

    return o_figure;
}

Figure creatDodecahedron(Color kleur){
    Figure d_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Figure I_figure = createIcosahedron(kleur);
    std::vector<Vector3D> I_points = I_figure.points;
    std::vector<Face> I_faces = I_figure.faces;

    for(int i = 0; i < 20; i++){
        Face face = I_faces[i];
        double x = 0;
        double y = 0;
        double z = 0;
        for(auto &index : face.point_indexes){
            x+= I_points[index].x;
            y+= I_points[index].y;
            z+= I_points[index].z;
        }
        x/=3;
        y/=3;
        z/=3;
        Vector3D center = Vector3D::point(x, y, z);
        points.push_back(center);

    }

    Face face1;
    face1.point_indexes.push_back(0);
    face1.point_indexes.push_back(1);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(3);
    face1.point_indexes.push_back(4);

    Face face2;
    face2.point_indexes.push_back(0);
    face2.point_indexes.push_back(5);
    face2.point_indexes.push_back(6);
    face2.point_indexes.push_back(7);
    face2.point_indexes.push_back(1);

    Face face3;
    face3.point_indexes.push_back(1);
    face3.point_indexes.push_back(7);
    face3.point_indexes.push_back(8);
    face3.point_indexes.push_back(9);
    face3.point_indexes.push_back(2);

    Face face4;
    face4.point_indexes.push_back(2);
    face4.point_indexes.push_back(9);
    face4.point_indexes.push_back(10);
    face4.point_indexes.push_back(11);
    face4.point_indexes.push_back(3);

    Face face5;
    face5.point_indexes.push_back(3);
    face5.point_indexes.push_back(11);
    face5.point_indexes.push_back(12);
    face5.point_indexes.push_back(13);
    face5.point_indexes.push_back(4);

    Face face6;
    face6.point_indexes.push_back(4);
    face6.point_indexes.push_back(13);
    face6.point_indexes.push_back(14);
    face6.point_indexes.push_back(5);
    face6.point_indexes.push_back(0);

    Face face7;
    face7.point_indexes.push_back(19);
    face7.point_indexes.push_back(18);
    face7.point_indexes.push_back(17);
    face7.point_indexes.push_back(16);
    face7.point_indexes.push_back(15);

    Face face8;
    face8.point_indexes.push_back(19);
    face8.point_indexes.push_back(14);
    face8.point_indexes.push_back(13);
    face8.point_indexes.push_back(12);
    face8.point_indexes.push_back(18);

    Face face9;
    face9.point_indexes.push_back(18);
    face9.point_indexes.push_back(12);
    face9.point_indexes.push_back(11);
    face9.point_indexes.push_back(10);
    face9.point_indexes.push_back(17);

    Face face10;
    face10.point_indexes.push_back(17);
    face10.point_indexes.push_back(10);
    face10.point_indexes.push_back(9);
    face10.point_indexes.push_back(8);
    face10.point_indexes.push_back(16);

    Face face11;
    face11.point_indexes.push_back(16);
    face11.point_indexes.push_back(8);
    face11.point_indexes.push_back(7);
    face11.point_indexes.push_back(6);
    face11.point_indexes.push_back(15);

    Face face12;
    face12.point_indexes.push_back(15);
    face12.point_indexes.push_back(6);
    face12.point_indexes.push_back(5);
    face12.point_indexes.push_back(14);
    face12.point_indexes.push_back(19);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);
    faces.push_back(face4);
    faces.push_back(face5);
    faces.push_back(face6);
    faces.push_back(face7);
    faces.push_back(face8);
    faces.push_back(face9);
    faces.push_back(face10);
    faces.push_back(face11);
    faces.push_back(face12);

    d_figure.points = points;
    d_figure.faces = faces;
    d_figure.color = kleur;

    return d_figure;



}

Figure createSphere(const int n, Color kleur){
    Figure s_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Figure i_figure = createIcosahedron(kleur);
    std::vector<Vector3D> i_points = i_figure.points;
    std::vector<Face> i_faces = i_figure.faces;

    std::vector<Face> new_faces;
    std::vector<Vector3D> new_points;
    if(n == 0){
        s_figure.points = i_points;
        s_figure.faces = i_faces;
        s_figure.color = kleur;
        return s_figure;
    }
    for(int i = 0; i < n; i++){
        std::vector<Face> temp_faces = i_faces;
        i_faces.clear();
        std::vector<Vector3D> temp_points = i_points;
        i_points.clear();
        for(auto &face : temp_faces){
            Vector3D p1 = temp_points[face.point_indexes[0]];
            Vector3D p2 = temp_points[face.point_indexes[1]];
            Vector3D p3 = temp_points[face.point_indexes[2]];

            Vector3D p4 = (p1 + p2)/2;
            Vector3D p5 = (p2 + p3)/2;
            Vector3D p6 = (p1 + p3)/2;


            int index1 = i_points.size();
            i_points.push_back(p1);
            int index2 = i_points.size();
            i_points.push_back(p2);
            int index3 = i_points.size();
            i_points.push_back(p3);
            int index4 = i_points.size();
            i_points.push_back(p4);
            int index5 = i_points.size();
            i_points.push_back(p5);
            int index6 = i_points.size();
            i_points.push_back(p6);

            Face new_face1;
            new_face1.point_indexes.push_back(index1);
            new_face1.point_indexes.push_back(index4);
            new_face1.point_indexes.push_back(index6);

            Face new_face2;
            new_face2.point_indexes.push_back(index4);
            new_face2.point_indexes.push_back(index2);
            new_face2.point_indexes.push_back(index5);

            Face new_face3;
            new_face3.point_indexes.push_back(index6);
            new_face3.point_indexes.push_back(index5);
            new_face3.point_indexes.push_back(index3);

            Face new_face4;
            new_face4.point_indexes.push_back(index4);
            new_face4.point_indexes.push_back(index5);
            new_face4.point_indexes.push_back(index6);

            i_faces.push_back(new_face1);
            i_faces.push_back(new_face2);
            i_faces.push_back(new_face3);
            i_faces.push_back(new_face4);


        }
        faces = i_faces;
        points = i_points;
    }

    for(auto &point : points){
        point.normalise();
    }
    s_figure.points = points;
    s_figure.faces = faces;
    s_figure.color = kleur;

    return s_figure;



}

Figure createCube(Color kleur){
    Figure c_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    Vector3D p1 = Vector3D::point(1, -1, -1);
    Vector3D p2 = Vector3D::point(-1, 1, -1);
    Vector3D p3 = Vector3D::point(1, 1, 1);
    Vector3D p4 = Vector3D::point(-1, -1, 1);
    Vector3D p5 = Vector3D::point(1, 1, -1);
    Vector3D p6 = Vector3D::point(-1, -1, -1);
    Vector3D p7 = Vector3D::point(1, -1, 1);
    Vector3D p8 = Vector3D::point(-1, 1, 1);

    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);
    points.push_back(p5);
    points.push_back(p6);
    points.push_back(p7);
    points.push_back(p8);

    Face face1;
    face1.point_indexes.push_back(0);
    face1.point_indexes.push_back(4);
    face1.point_indexes.push_back(2);
    face1.point_indexes.push_back(6);

    Face face2;
    face2.point_indexes.push_back(4);
    face2.point_indexes.push_back(1);
    face2.point_indexes.push_back(7);
    face2.point_indexes.push_back(2);

    Face face3;
    face3.point_indexes.push_back(1);
    face3.point_indexes.push_back(5);
    face3.point_indexes.push_back(3);
    face3.point_indexes.push_back(7);

    Face face4;
    face4.point_indexes.push_back(5);
    face4.point_indexes.push_back(0);
    face4.point_indexes.push_back(6);
    face4.point_indexes.push_back(3);

    Face face5;
    face5.point_indexes.push_back(6);
    face5.point_indexes.push_back(2);
    face5.point_indexes.push_back(7);
    face5.point_indexes.push_back(3);

    Face face6;
    face6.point_indexes.push_back(0);
    face6.point_indexes.push_back(5);
    face6.point_indexes.push_back(1);
    face6.point_indexes.push_back(4);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);
    faces.push_back(face4);
    faces.push_back(face5);
    faces.push_back(face6);

    c_figure.points = points;
    c_figure.faces = faces;
    c_figure.color = kleur;

    return c_figure;

}

Figure createCone(const int n, const double h, Color kleur){
    Figure c_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Vector3D top = Vector3D::point(0, 0, h);
    for(int i = 0; i < n; i++){
        double x = cos(2*PI*i/n);
        double y = sin(2*PI*i/n);
        Vector3D p = Vector3D::point(x, y, 0);
        points.push_back(p);
    }
    points.push_back(top);

    for(int i = 0; i <= n; i++){
        Face face;
        face.point_indexes.push_back(i);
        face.point_indexes.push_back((i+1)%n);
        face.point_indexes.push_back(n);
        faces.push_back(face);
    }

    Face face_n;
    for(int i = n-1; i>=0; i--){
        face_n.point_indexes.push_back(i);
    }
    faces.push_back(face_n);

    Face bottomFace;
    for(int i = 0; i < n; i++) {
        bottomFace.point_indexes.push_back(i);
    }
    faces.push_back(bottomFace);

    c_figure.points = points;
    c_figure.faces = faces;
    c_figure.color = kleur;

    return c_figure;


}

Figure createCylinder(const int n, const double h, Color kleur){
    Figure c_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    for(int i = 0; i < n; i++){
        double x = cos(2*3.14159265358979323846*i/n);
        double y = sin(2*3.14159265358979323846*i/n);
        Vector3D p1 = Vector3D::point(x, y, 0);
        Vector3D p2 = Vector3D::point(x, y, h);
        points.push_back(p1);
        points.push_back(p2);
    }

    for(int i = 0; i < n; i++){
        Face face1;
        face1.point_indexes.push_back(2*i);
        face1.point_indexes.push_back((2*i+2)%(2*n));
        face1.point_indexes.push_back((2*i+3)%(2*n));
        face1.point_indexes.push_back(2*i+1);
        faces.push_back(face1);

        Face face2;
        face2.point_indexes.push_back(2*i);
        face2.point_indexes.push_back(2*i+1);
        face2.point_indexes.push_back((2*i+3)%(2*n));
        face2.point_indexes.push_back((2*i+2)%(2*n));
        faces.push_back(face2);
    }

    Face topFace;
    for(int i = 0; i < 2*n; i++) {
        if(i % 2 == 1){
            topFace.point_indexes.push_back(i);
        }



    }
    faces.push_back(topFace);

    // Create bottom face
    Face bottomFace;
    for(int i = 0; i < 2*n; i++) {
        if(i % 2 == 0){
            bottomFace.point_indexes.push_back(i);
        }


    }
    faces.push_back(bottomFace);

    c_figure.points = points;
    c_figure.faces = faces;
    c_figure.color = kleur;

    return c_figure;

}

Figure createTorus(const int n, const int m, const double R, const double r, Color kleur){
    Figure t_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            Vector3D p = Vector3D::point((R+r*cos(2*j*PI/m))*cos(2*i*PI/n), (R+r*cos(2*j*PI/m))*sin(2*i*PI/n), r*sin(2*j*PI/m));
            points.push_back(p);
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            Face face1;
            face1.point_indexes.push_back(i*m+j);
            face1.point_indexes.push_back((i+1)%n*m+j);
            face1.point_indexes.push_back((i+1)%n*m+(j+1)%m);
            face1.point_indexes.push_back(i*m+(j+1)%m);
            faces.push_back(face1);
        }
    }

    t_figure.points = points;
    t_figure.faces = faces;
    t_figure.color = kleur;

    return t_figure;



}

Figure createBuckyBall(Color kleur) {
    Figure b_figure;
    std::vector<Vector3D> points;
    std::vector<Face> faces;

    // Creëer een icosahedron
    Figure I_figure = createIcosahedron(kleur);
    std::vector<Vector3D> I_points = I_figure.points;
    std::vector<Face> I_faces = I_figure.faces;

    // Voor elke driehoek van de icosahedron
    for (const auto& face : I_faces) {
        // Verdeel de driehoek in een gelijkzijdige zeshoek en drie driehoeken
        // Voeg de punten van de zeshoek toe aan de lijst van punten
        Vector3D center = (I_points[face.point_indexes[0]] + I_points[face.point_indexes[1]] + I_points[face.point_indexes[2]]) / 3.0;
        points.push_back(center);
        for (int i = 0; i < 3; ++i) {
            Vector3D v1 = I_points[face.point_indexes[i]];
            Vector3D v2 = I_points[face.point_indexes[(i + 1) % 3]];
            points.push_back((v1 + v2) / 2.0);
        }

        // Voeg de piramide met vijf zijden toe aan de lijst van gezichten
        Face pyramid_face;
        for (int i = 0; i < 4; ++i) {
            pyramid_face.point_indexes.push_back(points.size() - 4 + i);
        }
        faces.push_back(pyramid_face);
    }

    // Alleen het grondvlak van elke piramide behouden om de buckyball te construeren
    for (const auto& face : faces) {
        Face ground_face;
        ground_face.point_indexes.push_back(face.point_indexes[0]);
        ground_face.point_indexes.push_back(face.point_indexes[1]);
        ground_face.point_indexes.push_back(face.point_indexes[2]);
        faces.push_back(ground_face);
    }

    // Voeg de punten en gezichten toe aan de figuur
    b_figure.points = points;
    b_figure.faces = faces;
    b_figure.color = kleur;

    return b_figure;
}


void CalculateConstantes(Vector3D const& A, Vector3D const& B, Vector3D const& C, double& dzdx, double& dzdy, double& d){
    Vector3D vec_u = B - A;
    Vector3D vec_v = C - A;

    double w1 = vec_u.y*vec_v.z - vec_u.z*vec_v.y;
    double w2 = vec_u.z*vec_v.x - vec_u.x*vec_v.z;
    double w3 = vec_u.x*vec_v.y - vec_u.y*vec_v.x;

    double k = w1*A.x + w2*A.y + w3*A.z;

    dzdx = w1 / (-d*k);
    dzdy = w2 / (-d*k);
}


void draw_zbuf_triag(ZBuffer& buffer , img::EasyImage& image, Vector3D const& A, Vector3D const& B, Vector3D const& C, double& d, double& dx, double& dy, const img::Color& color){
    //Projected Points
    Point2D A_project{};
    A_project.x = (d*A.x/-A.z) + dx;
    A_project.y = (d*A.y/-A.z) + dy;

    Point2D B_project{};
    B_project.x = (d*B.x/-B.z) + dx;
    B_project.y = (d*B.y/-B.z) + dy;

    Point2D C_project{};
    C_project.x = (d*C.x/-C.z) + dx;
    C_project.y = (d*C.y/-C.z) + dy;

    int ymin = round(std::min({A_project.y, B_project.y, C_project.y}) + 0.5);
    int ymax = round(std::max({A_project.y, B_project.y, C_project.y}) - 0.5);

    vector<pair<Point2D, Point2D>> posible_snijpunten;
    //Emplace will make a pair and put it in the vector
    posible_snijpunten.emplace_back(A_project, B_project);
    posible_snijpunten.emplace_back(A_project, C_project);
    posible_snijpunten.emplace_back(B_project, C_project);


    double dzdx;
    double dzdy;
    CalculateConstantes(A,B,C, dzdx, dzdy, d);

    double xG = (A_project.x + B_project.x + C_project.x) / 3;
    double yG = (A_project.y + B_project.y + C_project.y) / 3;

    double zbw_G = (1/(3*A.z)) + (1/(3*B.z)) + (1/(3*C.z));

    for(int yi = ymin ; yi<= ymax; yi++){
        int x_L;
        int x_R;

        //Initialize variables
        double x_l_AB = std::numeric_limits<double>::infinity();
        double x_l_AC = std::numeric_limits<double>::infinity();
        double x_l_BC = std::numeric_limits<double>::infinity();

        double x_r_AB = -std::numeric_limits<double>::infinity();
        double x_r_AC = -std::numeric_limits<double>::infinity();
        double x_r_BC = -std::numeric_limits<double>::infinity();

        int count = 0;
        for(auto& lijnstuk : posible_snijpunten){
            double x_q = lijnstuk.second.x;
            double x_p = lijnstuk.first.x;
            double y_p = lijnstuk.first.y;
            double y_q = lijnstuk.second.y;

            if((yi-y_p)*(yi-y_q) <= 0 && y_p != y_q){
                double x_i =x_q + (x_p - x_q) * ((yi - y_q)/ (y_p - y_q));
                switch (count) {
                    case 0: {
                        x_l_AB = x_i;
                        x_r_AB = x_i;
                    }
                        break;
                    case 1: {
                        x_l_AC = x_i;
                        x_r_AC = x_i;
                    }
                        break;
                    case 2: {
                        x_l_BC = x_i;
                        x_r_BC = x_i;
                    }
                        break;

                }
                count++;

            }
        }
        x_L = round(min({x_l_AB,x_l_AC, x_l_BC}) + 0.5);
        x_R = round(max({x_r_AB, x_r_AC, x_r_BC}) - 0.5);



        for(int i = x_L; i<=x_R; i++){
            double zb_w = (1.0001 * zbw_G) + (i - xG)*dzdx + (yi - yG)*dzdy;
            if(zb_w < buffer.buffer[i][yi]){
                (image)(i, yi) = color;
                buffer.buffer[i][yi] = zb_w;

            }

        }



    }
}



void generateFractal(Figure& fig, Figures3D& fractal,const int nr_iterations, const double scale){

    if(nr_iterations == 0){
        return;
    }

    Matrix scala_matrix = scaleFigure(1/scale);

    for(int j = 0; j < nr_iterations; j++){
        Figures3D temp_fractal;
        for(auto& _fig : fractal){
            for(size_t i = 0; i < _fig.points.size(); i++){
                Figure fig_copy = _fig;
                applyTransformation(fig_copy, scala_matrix);

                Vector3D Pi = _fig.points[i];
                Vector3D Pi_acc = fig_copy.points[i];

                Vector3D trans_vector = Pi - Pi_acc;
                Matrix trans_matrix = translate(trans_vector);
                applyTransformation(fig_copy, trans_matrix);

                temp_fractal.push_back(fig_copy);
            }


        }
        fractal.clear();
        fractal.insert(fractal.end(), temp_fractal.begin(), temp_fractal.end());
        temp_fractal.clear();
    }



}

img::EasyImage generate_image(const ini::Configuration &configuration)
{

    if(configuration["General"]["type"].as_string_or_die() == "2DLSystem"){
        LParser::LSystem2D l_system;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());
        input_stream >> l_system;
        list<Line2D> line_s = drawLSystem(l_system, configuration);
        return draw2DLines(line_s, configuration["General"]["size"].as_int_or_die(), configuration);
    }



    else if(configuration["General"]["type"].as_string_or_die() == "Wireframe" || configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){
        Figures3D  figures;
        std::vector<double> eye_point = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D v_point = Vector3D::point(eye_point[0], eye_point[1], eye_point[2]);



        for(int i = 0; i <configuration["General"]["nrFigures"].as_int_or_die(); i++ ){
            Figure figuur{};
            Figures3D fractal;
            std::string fignum = "Figure"+ std::to_string(i);
            std::vector<double> clr = configuration[fignum]["color"].as_double_tuple_or_die();
            figuur.color.red = clr[0];
            figuur.color.green = clr[1];
            figuur.color.blue = clr[2];


            if(configuration[fignum]["type"].as_string_or_die() == "Tetrahedron"){
                figuur = createTetrahedron(figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Icosahedron"){
                figuur = createIcosahedron(figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Octahedron"){
                figuur = creatOctahedron(figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Dodecahedron"){
                figuur = creatDodecahedron(figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Sphere"){
                figuur = createSphere(configuration[fignum]["n"].as_int_or_die(), figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Cube"){
                figuur = createCube(figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Cone"){
                figuur = createCone(configuration[fignum]["n"].as_int_or_die(), configuration[fignum]["height"].as_double_or_die(), figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Cylinder"){
                figuur = createCylinder(configuration[fignum]["n"].as_int_or_die(), configuration[fignum]["height"].as_double_or_die(), figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "Torus"){
                figuur = createTorus(configuration[fignum]["n"].as_int_or_die(), configuration[fignum]["m"].as_int_or_die(), configuration[fignum]["R"].as_double_or_die(), configuration[fignum]["r"].as_double_or_die(), figuur.color);
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "BuckyBall"){
                figuur = createBuckyBall(figuur.color);
            }

            else if(configuration[fignum]["type"].as_string_or_die() == "3DLSystem"){
                LParser::LSystem3D l_system;
                std::ifstream input_stream(configuration[fignum]["inputfile"].as_string_or_die());
                input_stream >> l_system;

                figuur = draw3DLSystem(l_system, configuration, fignum);

            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalTetrahedron"){
                Figure fig = createTetrahedron(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());

            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalIcosahedron"){
                Figure fig = createIcosahedron(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalOctahedron"){
                Figure fig = creatOctahedron(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalDodecahedron"){
                Figure fig = creatDodecahedron(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalSphere"){
                Figure fig = createSphere(configuration[fignum]["n"].as_int_or_die(), figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalCube"){
                Figure fig = createCube(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }
            else if(configuration[fignum]["type"].as_string_or_die() == "FractalBuckyBall"){
                Figure fig = createBuckyBall(figuur.color);
                fractal.push_back(fig);
                generateFractal(fig, fractal, configuration[fignum]["nrIterations"].as_int_or_die(), configuration[fignum]["fractalScale"].as_double_or_die());
            }

            else if(configuration[fignum]["type"].as_string_or_die() == "LineDrawing"){
                for(int j =0; j< configuration[fignum]["nrLines"].as_int_or_die(); j++ ){
                    Face face{};
                    std::string linenum = "line"+ std::to_string(j);
                    std::vector<int> vec = configuration[fignum][linenum].as_int_tuple_or_die();

                    for(auto &num : vec){
                        face.point_indexes.push_back(num);
                    }
                    figuur.faces.push_back(face);
                }

                for(int k = 0; k < configuration[fignum]["nrPoints"].as_int_or_die(); k++){
                    std::string pointnum = "point"+ std::to_string(k);


                    std::vector<double> vec = configuration[fignum][pointnum].as_double_tuple_or_die();
                    Vector3D point = Vector3D::point(vec[0], vec[1], vec[2]);

                    figuur.points.push_back(point);

                }
            }
            else{
                continue;
            }



            Matrix scale_matrix = scaleFigure(configuration[fignum]["scale"].as_double_or_die());

            double degreeX = configuration[fignum]["rotateX"].as_double_or_die();
            Matrix rotateX_matrix = rotateX(degreesToRadians(degreeX));

            double degreeY = configuration[fignum]["rotateY"].as_double_or_die();
            Matrix rotateY_matrix = rotateY(degreesToRadians(degreeY));

            double degreeZ = configuration[fignum]["rotateZ"].as_double_or_die();
            Matrix rotateZ_matrix = rotateZ(degreesToRadians(degreeZ));

            std::vector<double> cen_ter = configuration[fignum]["center"].as_double_tuple_or_die();

            Vector3D center_vector = Vector3D::vector(cen_ter[0], cen_ter[1], cen_ter[2]);

            Matrix translate_matrix = translate(center_vector);

            Matrix eyePointTrans_matrix = eyePointTrans(v_point);

            Matrix total_matrix = scale_matrix * rotateX_matrix * rotateY_matrix * rotateZ_matrix * translate_matrix * eyePointTrans_matrix;
            applyTransformation(figuur, total_matrix);

            figures.push_back(figuur);

            for(auto& fig : fractal){
                applyTransformation(fig, total_matrix);
            }

            figures.insert(figures.end(), fractal.begin(), fractal.end());

        }

        list<Line2D> line_s = doProjection(figures);
        if(configuration["General"]["type"].as_string_or_die() == "ZBuffering"){
            double imageX;
            double imageY;
            double dx;
            double dy;
            double schaalfactor;

            Calculate(imageX, imageY, dx, dy, schaalfactor, line_s, configuration["General"]["size"].as_int_or_die(), configuration);

            img::EasyImage image(imageX, imageY);
            std::vector<double> bgcolor = configuration["General"]["backgroundcolor"];
            image.clear(img::Color(bgcolor[0]*255, bgcolor[1]*255, bgcolor[2]*255));
            ZBuffer buffer(imageX, imageY);


            for(auto& fig : figures){
                vector<Face> newFace;
                for(auto&face : fig.faces){
                    vector<Face> temp = triangulate(face);
                    newFace.insert(newFace.end(), temp.begin(), temp.end());
                }
                fig.faces = newFace;

            }


            for(auto& fig : figures){
                for(auto&face : fig.faces){
                    draw_zbuf_triag(buffer, image, fig.points[face.point_indexes[0]], fig.points[face.point_indexes[1]],fig.points[face.point_indexes[2]], schaalfactor, dx, dy, img::Color(fig.color.red*255, fig.color.green*255, fig.color.blue*255));
                }
            }
            return image;

        }
        else{
            return draw2DLines(line_s, configuration["General"]["size"].as_int_or_die(), configuration);
        };





    }


}




int main(int argc, char const* argv[])
{
    int retVal = 0;
    try
    {
        std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for(std::string fileName : args)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName);
                if (fin.peek() == std::istream::traits_type::eof()) {
                    std::cout << "Ini file appears empty. Does '" <<
                              fileName << "' exist?" << std::endl;
                    continue;
                }
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
