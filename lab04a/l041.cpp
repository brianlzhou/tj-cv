#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
#include <vector>
#include <stack>
#include <chrono>
using namespace std;

const int scale_row=400;
const int scale_column=400;
const int square_combinations=6;
const int square_pts=4;

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

RGB* pixel2D = new RGB[scale_row*scale_column];

void fill_pixel(RGB* imgcanvas, int x, int y) {
    if(x<0 || x>=400 || y<0 || y>=400) return;
    *(imgcanvas+y*scale_row+x)=RGB(0,0,0);
}

void fill_pixel_color(RGB* imgcanvas, int x, int y, int r, int g, int b) {
    if(x<0 || x>=400 || y<0 || y>=400) return;
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

        bool operator==(const Point& other) const {
            constexpr double epsilon = 1e-9;  // Precision threshold
            return std::fabs(x - other.x) < epsilon && std::fabs(y - other.y) < epsilon;
        }
};

class Square {
    private:
        unsigned long long r, c;
    
    public:
        Square(unsigned long long row, unsigned long long column) {
            r = row;
            c = column;
        }
    
        unsigned long long get_r() { return r; }
        unsigned long long get_c() { return c; }
    
        bool operator==(Square sqr) const {
            return r == sqr.get_r()   &&   c == sqr.get_c();
        }
};

class Pair {
    private:
        Point a=Point(0, 0);
        Point b=Point(0, 0);
        double dist;

    public:
        Pair() {
            a=Point(0, 0);
            b=Point(0, 0);
            dist=0.0;
        }

        Pair(Point pt1, Point pt2, double distance) {
            a=pt1;
            b=pt2;
            dist = distance;
        }

        Point get_pt1() { return a; }
        Point get_pt2() { return b; }
        double get_distance() { return dist; }
};

// Pair part1Pair, part2Pair, part3Pair, part4Pair;
// microseconds part1Time, part2Time, part3Time, part4Time;

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

// part 1 code
auto myPxls = new int[400][400];
Point p0;

double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p2.get_x()-p1.get_x(), 2)+pow(p2.get_y()-p1.get_y(), 2));
}

double crossProduct(const Point& O, const Point& A, const Point& B) {
    return (A.get_x() - O.get_x()) * (B.get_y() - O.get_y()) - (A.get_y() - O.get_y()) * (B.get_x() - O.get_x());
}

double polarAngle(const Point& reference, const Point& p) {
    return atan2(p.get_y() - reference.get_y(), p.get_x() - reference.get_x());
}

void findHull(vector<Point>& hull, const vector<Point>& points, const Point& pt1, const Point& pt2) {
    // find furtherst point
    int rightmostPoint=-1;
    double maxDistance=0.0;
    
    int pointssize=points.size();
    for (int i=0; i<pointssize; i++) {
        double distance=crossProduct(pt1, pt2, points[i]);
        if (distance>maxDistance) {
            rightmostPoint=i;
            maxDistance=distance;
        }
    }

    // add furtherst point to the convex hull
    if (rightmostPoint==-1) { return; }
    Point farthestPoint=points[rightmostPoint];
    hull.insert(hull.begin(), farthestPoint);

    // remove points in new triangle
    vector<Point> S1, S2;
    for (int i=0; i<pointssize; i++) {
        if (crossProduct(pt1, farthestPoint, points[i]) > 0)
            S1.push_back(points[i]);
        if (crossProduct(farthestPoint, pt2, points[i]) > 0)
            S2.push_back(points[i]);
    }

    // keep recurring
    findHull(hull, S1, pt1, farthestPoint);
    findHull(hull, S2, farthestPoint, pt2);
}

bool comparePointsYX(const Point& pt1, const Point& pt2) {
    return pt1.get_y()<pt2.get_y()||(pt1.get_y()==pt2.get_y()&&pt1.get_x()<pt2.get_x());
}

bool comparePointsPolarAngle(const Point& a, const Point& b, const Point& reference) {
    double angleA = polarAngle(reference, a);
    double angleB = polarAngle(reference, b);

    if (angleA == angleB) { return distance(reference, a) < distance(reference, b); } 
    else { return angleA < angleB; }
}

stack<Point> quickHull(const vector<Point>& points) {
    stack<Point> hullStack;
    int pointssize=points.size();
    if (pointssize<3) {
        for (int i=0; i<pointssize; i++) { hullStack.push(points[i]); }
        return hullStack;
    }

    // find the points w/ highest and lowest x coord
    int minX=0, maxX=0;
    for (int i=1; i<pointssize; i++) {
        if (points[i].get_x()<points[minX].get_x()) { minX=i; }
        if (points[i].get_x()>points[maxX].get_x()) { maxX=i; }
    }
    Point leftmost=points[minX], rightmost=points[maxX];
    vector<Point> leftSidePoints, rightSidePoints;
    for (int i=1; i<pointssize; i++) {
        if (crossProduct(leftmost, rightmost, points[i])>0) { leftSidePoints.push_back(points[i]); }
        else if (crossProduct(rightmost, leftmost, points[i])>0) { rightSidePoints.push_back(points[i]); }
    }

    // add these 2 points to your convex hull
    vector<Point> hullPoints;
    hullPoints.push_back(leftmost);
    hullPoints.push_back(rightmost);

    // recur on each side
    findHull(hullPoints, leftSidePoints, leftmost, rightmost);
    findHull(hullPoints, rightSidePoints, rightmost, leftmost);

    // add points to the convex hull
    int hullpointssize=hullPoints.size();
    for (int i=1; i<hullpointssize; i++) { hullStack.push(hullPoints[i]); }

    // sort
    Point referencePoint=*std::min_element(hullPoints.begin(), hullPoints.end(), comparePointsYX);
    sort(hullPoints.begin(), hullPoints.end(), [&referencePoint](const Point& a, const Point& b) { return comparePointsPolarAngle(a, b, referencePoint); });
    stack<Point> hull;
    for (const auto& p : hullPoints) { hull.push(p); }
    return hull;
}


int part1() {
    // generate 60 random points in the unit square (save them in a data structure of your chosing) and save them in points.txt in the same format as 3.4 
    // (each line has the x and y coordinate of a point separated by 2 spaces, use at least 20 digits precision)
    srand((unsigned int)time(NULL));
    string response;
    
    vector<Point> points;
    ofstream pointsfile;
    pointsfile.open("points.txt");
    pointsfile << setprecision(20);
    Point temp;
    
    int tempnumber=60;
    for (int i=0; i<tempnumber; i++) {
        temp = Point();
        pointsfile << temp.get_x() << "  " << temp.get_y() << endl;
        points.push_back(temp);
    }   
    pointsfile.close();
    
    int pointssize=points.size();
    if(pointssize != 60) { cerr << "The file does not contain 60 points"; }
    
    // for(const auto& point : points) {
    //     cout << setprecision(20);
    //     cout << "Point: (" << point.get_x() << ", " << point.get_y() << ")\n";
    // }

    // apply the QuickHull algorithm to find a convex hull that encloses all 60 points and the vertices are points from the set of 60 points you generated
    auto myStack = quickHull(points);

    // while (!myStack.empty()) {
    //     Point p = myStack.top();
    //     cout << "(" << p.get_x() << ", " << p.get_y() << ")" << std::endl;
    //     myStack.pop();
    // }

    // create a ppm of size 400x400 named quickhull.ppm in which you display a bold circle of radius 3 for each scaled point you generated 
    // (bold is obtained by doing a circle of radius 3, then radius 4), as well as connect the vertices you obtained for the convex hull 
    // (so in the image I should see 60 points and a convex hull with vertices circles of radius 3, use red color for the vertices, that are conected forming a convex hull)
    LineSegment l=LineSegment(0,0,0,0);
    Circle c=Circle(0,0,0);
    Circle c2=Circle(0,0,0);

    Point firstPoint=myStack.top();
    Point origPoint=myStack.top();
    myStack.pop();

    c=Circle(firstPoint.get_x()*scale_row, firstPoint.get_y()*scale_column, 3);
    c.draw_circle_color(pixel2D,255,0,0);
    c2=Circle(firstPoint.get_x()*scale_row, firstPoint.get_y()*scale_column, 4);
    c2.draw_circle_color(pixel2D,255,0,0);

    for (int i=0; i<pointssize; i++) {
        c=Circle(points[i].get_x()*scale_row, points[i].get_y()*scale_column, 3);
        c.draw_circle_color(pixel2D,0,0,0);
        c2=Circle(points[i].get_x()*scale_row, points[i].get_y()*scale_column, 4);
        c2.draw_circle_color(pixel2D,0,0,0);
    }

    c=Circle(firstPoint.get_x()*scale_row, firstPoint.get_y()*scale_column, 3);
    c.draw_circle_color(pixel2D,255,0,0);
    c2=Circle(firstPoint.get_x()*scale_row, firstPoint.get_y()*scale_column, 4);
    c2.draw_circle_color(pixel2D,255,0,0);

    while(!myStack.empty()) {
        Point curr=myStack.top();

        c=Circle(curr.get_x()*scale_row, curr.get_y()*scale_column, 3);
        c.draw_circle_color(pixel2D,255,0,0);
        c2=Circle(curr.get_x()*scale_row, curr.get_y()*scale_column, 4);
        c2.draw_circle_color(pixel2D,255,0,0);

        l=LineSegment(scale_row*firstPoint.get_x(), scale_column*firstPoint.get_y(), scale_row*curr.get_x(), scale_column*curr.get_y());
        l.draw_line(pixel2D);

        firstPoint=curr;
        myStack.pop();
    }
    l=LineSegment(scale_row*firstPoint.get_x(), scale_column*firstPoint.get_y(), scale_row*origPoint.get_x(), scale_column*origPoint.get_y());
    l.draw_line(pixel2D);
    
    ofstream outfile;
    srand(0);
    outfile.open("quickhull.ppm");
    outfile << "P3\n" << 400 << " " << "400\n" << "255\n";
    RGB temprgb;
    for (int r=0; r<scale_row; r++) {
        for (int c=0; c<scale_column; c++) {
            temprgb=pixel2D[r*scale_column+c];
            outfile << temprgb.get_red() << " " << temprgb.get_green() << " " << temprgb.get_blue() << " ";
        }
        outfile << "\n";
    }
    outfile.close();
    return 0;
}

int main() {
    cout << "test";
    part1();
}