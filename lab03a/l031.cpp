#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
using namespace std;

const int scale_row=800;
const int scale_column=800;
const int square_combinations=6;
const int square_pts=4;

// auto pixel2D = new bool[800][800];
auto pixel2D = new bool[800][800];

class RGB {
    private:
        int rval=255, gval = 255, bval = 255;

    public:
        RGB() {
            rval=255;
            gval=255;
            bval=255;
        }

        RGB(int rv, int gv, int bv) {
            rval=rv;
            gval=gv;
            bval=bv;
        }

        int get_red() { return rval; }
        int get_green() { return gval; }
        int get_blue() { return bval; }
};

void fill_pixel(RGB* imgcanvas, int x, int y) {
    if(x<0 || x>=800 || y<0 || y>=800) return;
    *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
}

void fill_pixel_color(RGB* imgcanvas, int x, int y, int r, int g, int b) {
    if(x<0 || x>=800 || y<0 || y>=800) return;
    *(imgcanvas+y*scale_row+x)=RGB(r,g,b);
}

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
    
        double get_x() const { return x; }
        double get_y() const { return y; }
};

class LineSegment {
    private:
        int x1, y1, x2, y2;

    public:
        LineSegment(int xA, int yA, int xB, int yB) {
            x1 = xA;
            y1 = yA;
            x2 = xB;
            y2 = yB;
        }

        void draw_line_XY(RGB* imgcanvas) {
            if (x1==x2) {
                if (y2>y1) {
                    for (int idx=y1; idx<=y2; idx++) {
                        *(imgcanvas+idx*scale_row+x1) = RGB(0,0,0);
                    }
                } else {
                    for (int idx=y2; idx<=y1; idx++) {
                        *(imgcanvas+idx*scale_row+x1) = RGB(0,0,0);
                    }
                }
            } else {
                if (x2>x1) {
                    for (int idx=x1; idx<=x2; idx++) {
                        *(imgcanvas+y1*scale_row+idx) = RGB(0,0,0);
                    }
                } else {
                    for (int idx=x2; idx<=x1; idx++) {
                        *(imgcanvas+y1*scale_row+idx) = RGB(0,0,0);
                    }
                }
            }
        }

        void draw_line_xP(int x1, int y1, int x2, int y2, RGB* imgcanvas){
            int dX=x2-x1,
                dY=y2-y1,
                y=y1,
                dYdXdiff=dY-dX;
            for (int x=x1; x<=x2-1; x++) {
                *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
                if (dYdXdiff>0) {
                    y += 1;
                    dYdXdiff -= dX;
                }
                dYdXdiff+=dY;
            }
//             int mx = 0;
        }

        void draw_line_xN(int x1, int y1, int x2, int y2, RGB* imgcanvas) {
            int dX=x2-x1,
                dY=y2-y1,
                y=y1,
                dYdXdiff=dX-dY;
            for (int x=x1; x<=x2-1; x++) {
                *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
                if (dYdXdiff>0) {
                    y-=1;
                    dYdXdiff-=dX;
                }
                dYdXdiff+=abs(dY);
            }
        }

        void draw_line_yP(int x1, int y1, int x2, int y2, RGB* imgcanvas) {
            int dX=x2-x1, 
                dY=y2-y1, 
                x=x1, 
                dYdXdiff=dX-dY;
            for (int y=y1; y<=y2-1; y++) {
                *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
                if (dYdXdiff > 0) {
                    x+=1;
                    dYdXdiff-=dY;
                }
                dYdXdiff+=dX;
            }
        }

        void draw_line_yN(int x1, int y1, int x2, int y2, RGB* imgcanvas) {
            int dX=x2-x1, 
                dY=y2-y1, 
                x=x1, 
                dYdXdiff=dX-dY;
            for (int y=y1; y<=y2-1; y++) {
                *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
                if (dYdXdiff>0) {
                    x-=1;
                    dYdXdiff-=dY;
                }
                dYdXdiff+=abs(dX);
            }
        }

        void draw_line(RGB* imgcanvas) {
            if (x1==x2||y1==y2) {
                draw_line_XY(imgcanvas);
            } else {
                int dx=x2-x1;
                int dy=y2-y1;

                if(dx>=dy)
                    if(dy<0 && dx<0) draw_line_yP(x2,y2,x1,y1,imgcanvas);
                    else if(dy < 0) 
                        if(abs(dx) > abs(dy))
                            draw_line_xN(x1,y1,x2,y2,imgcanvas);
                        else
                            draw_line_yN(x2,y2,x1,y1,imgcanvas);
                    else draw_line_xP(x1,y1,x2,y2,imgcanvas);
                else if(dy<0 && dx<0) draw_line_xP(x2,y2,x1,y1,imgcanvas);
                else if(dx < 0) 
                    if(abs(dx) < abs(dy))
                        draw_line_yN(x1,y1,x2,y2,imgcanvas);
                    else
                        draw_line_xN(x2,y2,x1,y1,imgcanvas);
                else draw_line_yP(x1,y1,x2,y2,imgcanvas);
            }
        }
};

class Circle {
    private:
        int centerX=0;
        int centerY=0;
        int radius=10;

    public:
        Circle(int x, int y, int r) {
            centerX=x;
            centerY=y;
            radius=r;
        }

        void draw_circle(RGB* imgcanvas) {
            int x, y, xmax, y2, y2_new, ty;
            
            xmax = (int) (radius * 0.70710678+0.5);

            y=radius;
            y2 = y * y;
            ty = (2 * y) - 1;
            y2_new = y2;

            for (x = 0; x <= xmax + 1; x++) {
                if ((y2 - y2_new) >= ty) {
                    y2 -= ty;
                    y -= 1;
                    ty -= 2;
                }
                fill_pixel(imgcanvas, centerX+x, centerY-y);
                fill_pixel(imgcanvas, centerX-x, centerY-y);
                fill_pixel(imgcanvas, centerX+x, centerY+y);
                fill_pixel(imgcanvas, centerX-x, centerY+y);
                fill_pixel(imgcanvas, centerX+y, centerY-x);
                fill_pixel(imgcanvas, centerX-y, centerY-x);
                fill_pixel(imgcanvas, centerX+y, centerY+x);
                fill_pixel(imgcanvas, centerX-y, centerY+x);
                y2_new -= (2 * x) - 3;
            }
        }

        void draw_circle_color(RGB* imgcanvas, int r, int g, int b) {
            int x, y, xmax, y2, y2_new, ty;
            
            xmax = (int) (radius * 0.70710678+0.5);

            y=radius;
            y2 = y * y;
            ty = (2 * y) - 1;
            y2_new = y2;

            for (x = 0; x <= xmax + 1; x++) {
                if ((y2 - y2_new) >= ty) {
                    y2 -= ty;
                    y -= 1;
                    ty -= 2;
                }
                fill_pixel_color(imgcanvas, centerX+x, centerY-y, r, g, b);
                fill_pixel_color(imgcanvas, centerX-x, centerY-y, r, g, b);
                fill_pixel_color(imgcanvas, centerX+x, centerY+y, r, g, b);
                fill_pixel_color(imgcanvas, centerX-x, centerY+y, r, g, b);
                fill_pixel_color(imgcanvas, centerX+y, centerY-x, r, g, b);
                fill_pixel_color(imgcanvas, centerX-y, centerY-x, r, g, b);
                fill_pixel_color(imgcanvas, centerX+y, centerY+x, r, g, b);
                fill_pixel_color(imgcanvas, centerX-y, centerY+x, r, g, b);
                y2_new -= (2 * x) - 3;
            }
        }
};

// part 0 code
void part0() {
//     g++ -std=c++11 -Wall -o hello hello.cpp
    srand((unsigned int)time(NULL));
    string response;
    
    cout << "If you would like to generate 60 random points, please type 'yes': ";
    cin >> response;
    
    if(response == "yes") {
        ofstream pointsfile;
        pointsfile.open("points.txt");
        pointsfile << setprecision(17);
        Point temp;
        for (int i=0; i<60; i++) {
            temp = Point();
            pointsfile << temp.get_x() << " " << temp.get_y() << endl;
        }   
        pointsfile.close();
    } else {
        std::cout << "You did not type 'yes'." << std::endl;
    }
}

// part 1 code
double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p2.get_x()-p1.get_x(), 2)+pow(p2.get_y()-p1.get_y(), 2));
}

int part1() {
    // setup points ppm file
    ofstream outfile;
    outfile.open("points.ppm");
    outfile << "P3 " << scale_column << " " << scale_row << " " << "255" << endl;

    RGB* pixel2D=new RGB[scale_column*scale_row];
    for (int c=0; c<scale_column; c++) {
        for (int r=0; r<scale_row; r++) {
            *(pixel2D+r*scale_column+c)=RGB();
        }
    }

    // Part a: Read the points that are in the points.txt file (file either generated by part0 or the file is already there) and save them in a list
    list<Point> points;
    ifstream file("points.txt");
    if (file.fail()) return -1;
    
    // read the points
    double x, y;
    while (file >> x >> y) {
        points.push_back(Point(x, y));
    }
    file.close();
    
    // check if something went wrong
    if(points.size() != 60) { cerr << "The file does not contain 60 points"; }
    
    // display the points
    // for(const auto& point : points) {
    //     cout << setprecision(17);
    //     cout << "Point: (" << point.get_x() << ", " << point.get_y() << ")\n";
    // }
    
    // Part b: solve the closest-pair problem in the unit square using the"brute force" approach discussed in class. (use at least one iterator over the list)
    double min_distance=10000000000000000.1;
    Point min_pt1=Point(0,0), min_pt2=Point(0,0);

    for (auto pt1=points.cbegin(); pt1!=points.cend(); ++pt1) {
        for (auto pt2=next(pt1); pt2!=points.cend(); ++pt2) {
            double current_distance=distance(*pt1, *pt2);
            if (current_distance<min_distance) {
                min_distance=current_distance;
                min_pt1=*pt1;
                min_pt2=*pt2;
            }
        }
    }

    // Part c: create points.ppm file in which you draw a bold black circle of radius 3 for every point you generated and you draw a bold red circle of radius 3 for the 2 closest points you found (you may make the circles bold by drawig a circle of radius 2 and a circle of radius 3) 
    for (auto pt=points.cbegin(); pt!=points.cend(); ++pt) {
        if ((pt->get_x()==min_pt1.get_x() && pt->get_y()==min_pt1.get_y()) || (pt->get_x()==min_pt2.get_x() && pt->get_y()==min_pt2.get_y())) {
            Circle c(pt->get_x()*scale_row, pt->get_y()*scale_column, 3);
            c.draw_circle_color(pixel2D,255,0,0);
            Circle c2(pt->get_x()*scale_row, pt->get_y()*scale_column, 2);
            c2.draw_circle_color(pixel2D,255,0,0);
        } else{
            Circle c(pt->get_x()*scale_row, pt->get_y()*scale_column, 3);
            c.draw_circle(pixel2D);
            Circle c2(pt->get_x()*scale_row, pt->get_y()*scale_column, 2);
            c2.draw_circle(pixel2D);
        }
    }

    RGB temprgb;
    for (int r=0; r<scale_row; r++) {
        for (int c=0; c<scale_column; c++) {
            temprgb=pixel2D[r*scale_column+c];
            outfile << temprgb.get_red() << " " << temprgb.get_green() << " " << temprgb.get_blue() << " ";
        }
        outfile << "\n";
    }
    outfile.close();

    // Part d: display on the screen the 2 closest points and the distance between them (so the minimum distance)
    cout << setprecision(17);
    cout << "The closest pair of points are: \n";
    cout << "Point 1: (" << min_pt1.get_x() << ", " << min_pt1.get_y() << ")\n";
    cout << "Point 2: (" << min_pt2.get_x() << ", " << min_pt2.get_y() << ")\n";
    cout << "The distance between them is: " << min_distance << "\n";

    return 0;
}

int main() {
    part0();
    part1();
}