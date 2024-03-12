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
        if(i== '+' || i == '-' || i == '(' || i == ')'){
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

    for(int i = 0; i < l_system.get_nr_iterations(); i++){
        std::string newstring = Replace(cur_string, l_system);
        cur_string = newstring;
    }

    return cur_string;
}

double degreesToRadians(double& degrees) {
    double PI = 3.14159265358979323846;
    return degrees * (PI / 180.0);
}

img::EasyImage draw2DLines(const Lines2D &lines,const int size, const ini::Configuration &configuration){
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


    double imageX = size * (xrange/ fmax(xrange,yrange));
    double imageY = size * (yrange/ fmax(xrange,yrange));


    double schaalfactor = 0.95 * (imageX/xrange);


    Lines2D newLines = lines;

    mult_cords_with_SF(schaalfactor, newLines);

    double DCx = schaalfactor * ((xmin + xmax)/2);
    double DCy = schaalfactor * ((ymin + ymax)/2);

    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;



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

    for(auto &line : newLines){

        image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y,img::Color((line.color.red)*255,(line.color.green)*255, (line.color.blue)*255));

    }
    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system){
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

img::EasyImage generate_image(const ini::Configuration &configuration)
{

    if(configuration["General"]["type"].as_string_or_die() == "2DLSystem"){
        LParser::LSystem2D l_system;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());
        input_stream >> l_system;
        return draw2DLines(drawLSystem(l_system), configuration["General"]["size"].as_int_or_die(), configuration);
    }



    else if(configuration["General"]["type"].as_string_or_die() == "Wireframe"){
        Figures3D  figures;
        std::vector<double> eye_point = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D v_point = Vector3D::point(eye_point[0], eye_point[1], eye_point[2]);



        for(int i = 0; i <configuration["General"]["nrFigures"].as_int_or_die(); i++ ){
            Figure figuur{};
            std::string fignum = "Figure"+ std::to_string(i);
            std::vector<double> clr = configuration[fignum]["color"].as_double_tuple_or_die();
            figuur.color.red = clr[0];
            figuur.color.green = clr[1];
            figuur.color.blue = clr[2];


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


            Matrix scale_matrix = scaleFigure(configuration[fignum]["scale"].as_double_or_die());
            applyTransformation(figuur, scale_matrix);


            double degreeX = configuration[fignum]["rotateX"].as_double_or_die();
            Matrix rotateX_matrix = rotateX(degreesToRadians(degreeX));
            applyTransformation(figuur, rotateX_matrix);

            double degreeY = configuration[fignum]["rotateY"].as_double_or_die();
            Matrix rotateY_matrix = rotateY(degreesToRadians(degreeY));
            applyTransformation(figuur, rotateY_matrix);

            double degreeZ = configuration[fignum]["rotateZ"].as_double_or_die();
            Matrix rotateZ_matrix = rotateZ(degreesToRadians(degreeZ));
            applyTransformation(figuur, rotateZ_matrix);

            std::vector<double> cen_ter = configuration[fignum]["center"].as_double_tuple_or_die();

            Vector3D center_vector = Vector3D::vector(cen_ter[0], cen_ter[1], cen_ter[2]);

            Matrix translate_matrix = translate(center_vector);
            applyTransformation(figuur, translate_matrix);


            Matrix eyePointTrans_matrix = eyePointTrans(v_point);
            applyTransformation(figuur, eyePointTrans_matrix);


            figures.push_back(figuur);

        }

        return draw2DLines(doProjection(figures), configuration["General"]["size"].as_int_or_die(), configuration);


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
