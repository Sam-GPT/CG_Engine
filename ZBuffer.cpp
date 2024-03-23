//
// Created by USER on 20/03/2024.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int w_idth,const int h_eight){

    width = w_idth;
    height = h_eight;

    buffer.resize(width, std::vector<double>(height, std::numeric_limits<double>::infinity()));




}

ZBuffer::~ZBuffer() {
    buffer.clear();
}
