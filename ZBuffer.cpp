//
// Created by USER on 20/03/2024.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width,const int height){

    if (width <= 0 || height <= 0 || width > 10000 || height > 10000) {
        cout<<"Width : "<<width<<" Height : "<<height<<endl;
        throw std::invalid_argument("Width and height must be positive and not excessively large.");
    }
    this->width = width;
    this->height = height;

    buffer.resize(width, std::vector<double>(height, std::numeric_limits<double>::infinity()));
    for(int i = 0; i<width; i++){
        for(int j = 0; j < height; j++){
            buffer[i][j] = std::numeric_limits<double>::infinity();
        }
    }




}

ZBuffer::~ZBuffer() {
    buffer.clear();
}
