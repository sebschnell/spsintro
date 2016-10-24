#ifndef POLYGON_H
#define POLYGON_H
#include <vector>
#include "utilities.h"

class Polygon {
  public:
    // constructors
    Polygon() = default;
    Polygon(const std::vector<double> &rx, const std::vector<double> &ry) : x(rx), y(ry) {}
    // function members
    std::vector<double> get_x() const { return x; }
    std::vector<double> get_y() const { return y; }
    virtual double area() const { return calc_poly_area(x, y); }
    virtual ~Polygon() = default; // dynamic binding for the destructor
  private:
    // data members
    std::vector<double> x;
    std::vector<double> y;
};

class Rectangle : Polygon {
  public:
    Rectangle() = default;
    Rectangle(const std::vector<double> &rx, const std::vector<double> &ry,
              const double &w, const double &h) :
      Polygon(rx, ry), width(w), height(h) {}
    double get_w() const {return width;}
    double get_h() const {return height;}
    double area() const override {return width*height;}
  private:
    double width = 0.0;
    double height = 0.0;
};

class Triangle : Polygon {
  public:
    Triangle() = default;
    Triangle(const std::vector<double> &rx, const std::vector<double> &ry,
              double b, double h) :
      Polygon(rx, ry), base(w), height(h) {}
    double get_b() const {return base;}
    double get_h() const {return height;}
    double area() const override {return base*height/2;};
  private:
    double base = 0.0;
    double height = 0.0;
};
#endif
