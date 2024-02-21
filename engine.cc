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

using Lines2D = std::list<Line2D>;
void getXmin(Line2D& line, double& cur_min){
    if(line.p1.x < cur_min){
        cur_min = line.p1.x;
    }
    else if(line.p2.x < cur_min){
        cur_min = line.p2.x;
    }
}

void getXmax(Line2D& line, double& cur_max){
    if(line.p1.x > cur_max){
        cur_max = line.p1.x;
    }
    else if(line.p2.x > cur_max){
        cur_max= line.p2.x;
    }
}

void getYmin(Line2D& line, double& cur_min){
    if(line.p1.y < cur_min){
        cur_min = line.p1.y;
    }
    else if(line.p2.y < cur_min){
        cur_min = line.p2.y;
    }
}

void getYmax(Line2D& line, double& cur_max){
    if(line.p1.y > cur_max){
        cur_max = line.p1.y;
    }
    else if(line.p2.y > cur_max){
        cur_max= line.p2.y;
    }
}
Lines2D mult_cords_with_SF(double& d, Lines2D lines){
    for(auto line : lines){
        line.p1.x *=d;
        line.p1.y *= d;
        line.p2.x *=d;
        line.p2.y *= d;
    }
    return lines;
}
img::EasyImage draw2DLines(const Lines2D &lines,const int size){
    double xmin = lines.begin()->p1.x;
    double xmax = lines.begin()->p1.x;
    double ymin = lines.begin()->p1.y;
    double ymax = lines.begin()->p1.y;
    for(auto line : lines){
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


    Lines2D newLines = (schaalfactor, lines);

    double DCx = schaalfactor * ((xmin + xmax)/2);
    double DCy = schaalfactor * ((ymin + ymax)/2);

    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;

    for(auto &line : newLines){
        line.p1.x+=dx;
        line.p1.x= round(line.p1.x);
        std::cout <<line.p1.x<<std::endl;
        line.p2.x+=dx;
        line.p2.x=round(line.p2.x);
        std::cout <<line.p2.x<<std::endl;

        line.p1.y+=dy;
        line.p1.y=round(line.p1.y);
        std::cout <<line.p1.y<<std::endl;
        line.p2.y+=dy;
        line.p2.y= round(line.p2.y);
        std::cout <<line.p2.y<<std::endl;
    }

    img::EasyImage image(imageX, imageY);

    for(auto &line : newLines){
        image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y,img::Color(line.color.red*255, line.color.blue*255, line.color.green*255));

    }

    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system){
    Lines2D lines;
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    Point2D p1;
    p1.x = -100;
    p1.y = -110;

    Point2D p2;
    p2.x = 20;
    p2.y = 50;

    Line2D l1;
    l1.p1 = p1;
    l1.p2 = p2;
    l1.color.red = 0.2;
    l1.color.green = 0.1;
    l1.color.blue = 0.5;

    Lines2D lines;
    lines.push_back(l1);
    return draw2DLines(lines, 500);
//    return img::EasyImage(100,100);
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
