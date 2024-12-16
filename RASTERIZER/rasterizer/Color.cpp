#include <iomanip>
#include "Color.h"

using namespace std;

Color::Color() {
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

std::ostream &operator<<(std::ostream &os, const Color &c)
{
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}

Color Color::operator+(const Color &color1)
{
    Color color;
    color.r = this->r + color1.r;
    color.g = this->g + color1.g;
    color.b = this->b + color1.b;
    return color;
}

Color Color::operator-(const Color &color1)
{
    Color color;
    color.r = this->r - color1.r;
    color.g = this->g - color1.g;
    color.b = this->b - color1.b;
    return color;
}

Color Color::operator*(double variable)
{
    Color color;
    color.r = this->r * variable;
    color.g = this->g * variable;
    color.b = this->b * variable;
    return color;
}

Color Color::operator/(double variable)
{
    Color color;
    color.r = this->r / variable;
    color.g = this->g / variable;
    color.b = this->b / variable;
    return color;
}

Color Color::operator+=(const Color &color1)
{
    this->r += color1.r;
    this->g += color1.g;
    this->b += color1.b;
    return *this;
}

Color Color::round()
{
    Color color;
    color.r = int(this->r + 0.5);
    color.g = int(this->g + 0.5);
    color.b = int(this->b + 0.5);
    return color;
}