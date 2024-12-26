#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
using namespace std;

class Point {
    private:
        double x=0.0;
        double y=0.0;

    public:
        Point() {
            x = double(rand()) / RAND_MAX;
            y = double(rand()) / RAND_MAX;
        }

        Point(double x_input, double y_input) {
            x = x_input;
            y = y_input;
        }
    
        double get_x() { return x; }
        double get_y() { return y; }
};

double get_triangle_area(Point a, Point b, Point c) {
    double xA=a.get_x(), yA=a.get_y();
    double xB=b.get_x(), yB=b.get_y();
    double xC=c.get_x(), yC=c.get_y();
    double area_abc=0.5*(abs(xA*(yB-yC)+xB*(yC-yA)+xC*(yA-yB)));
    return area_abc;
};

bool calc_inside_point(Point a, Point b, Point c, Point newpoint) {
    double area_ABC = get_triangle_area(a, b, c);
    double area_ABD = get_triangle_area(a, b, newpoint);
    double area_BCD = get_triangle_area(b, c, newpoint);
    double area_ACD = get_triangle_area(a, c, newpoint);
    if (abs(area_ABC-(area_ABD+area_BCD+area_ACD))<0.0001) {
        return true;
    }
    return false;
};

void part1() {
    ofstream outfile;
    outfile.open("log.txt");

    srand(time(0));

    Point a = Point();
    double xA = a.get_x(), yA = a.get_y();
    Point b = Point();
    double xB = b.get_x(), yB = b.get_y();
    Point c = Point();
    double xC = c.get_x(), yC = c.get_y();
    
    //     replace duplicate points
    while ((xA==xB) && (yA==yB)) {
        b = Point();
        xB = b.get_x(), yB = b.get_y();
    }
    while ((xC==xA && yC==yA) || (xC==xB && yC==yB)) {
        c = Point();
        xC = c.get_x(), yC = c.get_y();
    }
    
    //     create a new point
    Point newpoint = Point();
    double xNewpoint = newpoint.get_x(), yNewpoint = newpoint.get_y();
    while ((xNewpoint==xA && yNewpoint==yA) || (xNewpoint==xB && yNewpoint==yB) || (xNewpoint==xC && yNewpoint==yC)) {
        newpoint = Point();
        xNewpoint = newpoint.get_x(), yNewpoint = newpoint.get_y();
    }
    outfile << fixed << setprecision(17);
    outfile << "Point A: (" << xA << "," << yA << ")" << endl;
    outfile << "Point B: (" << xB << "," << yB << ")" << endl;
    outfile << "Point C: (" << xC << "," << yC << ")" << endl;
    
    //     find a valid new point
    outfile << "testing point: (" << xNewpoint << "," << yNewpoint << ") \n";
    while (calc_inside_point(a,b,c,newpoint) || calc_inside_point(a,b,newpoint,c) || calc_inside_point(a,newpoint,c,b) || calc_inside_point(newpoint,c,b,a)) {
        newpoint = Point();
        xNewpoint = newpoint.get_x(), yNewpoint = newpoint.get_y();
        outfile << "testing point (" << xNewpoint << "," << yNewpoint << ") \n";
    }
    outfile.close();
    
    //     write valid points.txt file
    ofstream outfile2;
    outfile2.open("points.txt");

    outfile2 << fixed << setprecision(17);
    outfile2 << "(" << xA << "," << yA << ") , ";
    outfile2 << "(" << xB << "," << yB << ") , ";
    outfile2 << "(" << xC << "," << xC << ") , ";
    outfile2 << "(" << xNewpoint << "," << yNewpoint << ")";
    outfile2.close();
}

int main() {
    part1();
}