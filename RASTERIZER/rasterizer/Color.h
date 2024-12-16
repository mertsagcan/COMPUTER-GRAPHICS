#ifndef __COLOR_H__
#define __COLOR_H__
#include <iostream>
#include <cmath>

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    friend std::ostream &operator<<(std::ostream &os, const Color &c);

    Color round();
    
    //Operator funcions. (+, -, *, /))
    Color operator+(const Color &color1);
    Color operator-(const Color &color1);
    Color operator*(double variable);
    Color operator/(double variable);
    Color operator+=(const Color &color1);
    

};

#endif