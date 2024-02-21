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

std::string Replace(std::string string, const LParser::LSystem2D &l_system){
    std::string newstring = "";

    for(auto i : string){
        if(i== '+' || i == '-' || i == '[' || i == ']'){
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
    /*
    std::ofstream outFile("output.txt");
    outFile <<cur_string <<std::endl;
    outFile.close();

    */
    return cur_string;
}

double degreesToRadians(double& degrees) {
    double PI = 3.14159265358979323846;
    return degrees * (PI / 180.0);
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

    std::cout << imageX <<" " << imageY<<std::endl;
    double schaalfactor = 0.95 * (imageX/xrange);


    Lines2D newLines = lines;

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
    for(auto &line : newLines){

        image.draw_line(line.p1.x, line.p1.y, line.p2.x, line.p2.y,img::Color(line.color.red*255, line.color.blue*255, line.color.green*255));

    }
    return image;
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system){
    Lines2D lines;

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
            newLine.color.red = 0.2;
            newLine.color.green = 0.1;
            newLine.color.blue = 0.5;

            lines.push_back(newLine);
            cur_position.x = new_position.x;
            cur_position.y = new_position.y;


        }
        else if(i == '+'){
            cur_angle+= angle;
        }
        else if(i == '-'){
            cur_angle-=angle;
        }
        else{
            continue;
        }
    }


    return lines;



}

img::EasyImage generate_image(const ini::Configuration &configuration)
{

    LParser::LSystem2D l_system;
    std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());
    input_stream >> l_system;



    return draw2DLines(drawLSystem(l_system), configuration["General"]["size"].as_int_or_die());
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
