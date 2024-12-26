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
#include <chrono>
#include <unordered_map> 
using namespace std;
using namespace std::chrono;

const int scale_row=800;
const int scale_column=800;
const int square_combinations=6;
const int square_pts=4;

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

Pair part1Pair, part2Pair, part3Pair, part4Pair;
microseconds part1Time, part2Time, part3Time, part4Time;

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
    
    int tempnumber=2000000;
    
    cout << "If you would like to generate " << tempnumber << " random points, please type 'yes': ";
    cin >> response;
    
    if(response == "yes") {
        ofstream pointsfile;
        pointsfile.open("points.txt");
        pointsfile << setprecision(17);
        Point temp;
        for (int i=0; i<tempnumber; i++) {
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

// int part1() {
//     // setup points ppm file
//     ofstream outfile;
//     outfile.open("points.ppm");
//     outfile << "P3 " << scale_column << " " << scale_row << " " << "255" << endl;

//     RGB* pixel2D=new RGB[scale_column*scale_row];
//     for (int c=0; c<scale_column; c++) {
//         for (int r=0; r<scale_row; r++) {
//             *(pixel2D+r*scale_column+c)=RGB();
//         }
//     }

//     // Part a: Read the points that are in the points.txt file (file either generated by part0 or the file is already there) and save them in a list
//     list<Point> points;
//     ifstream file("points.txt");
//     if (file.fail()) return -1;
    
//     // read the points
//     double x, y;
//     while (file >> x >> y) {
//         points.push_back(Point(x, y));
//     }
//     file.close();
    
//     // check if something went wrong
//     // if(points.size() != 60) { cerr << "The file does not contain 60 points"; }
    
//     // display the points
//     // for(const auto& point : points) {
//     //     cout << setprecision(17);
//     //     cout << "Point: (" << point.get_x() << ", " << point.get_y() << ")\n";
//     // }
    
//     // Part b: solve the closest-pair problem in the unit square using the"brute force" approach discussed in class. (use at least one iterator over the list)
//     auto startTime=high_resolution_clock::now();
//     double min_distance=10000000000000000.1;
//     Point min_pt1=Point(0,0), min_pt2=Point(0,0);

//     for (auto pt1=points.cbegin(); pt1!=points.cend(); ++pt1) {
//         for (auto pt2=next(pt1); pt2!=points.cend(); ++pt2) {
//             double current_distance=distance(*pt1, *pt2);
//             if (current_distance<min_distance) {
//                 min_distance=current_distance;
//                 min_pt1=*pt1;
//                 min_pt2=*pt2;
//             }
//         }
//     }

//     auto stopTime=high_resolution_clock::now();
//     auto elapsedTime=duration_cast<microseconds>(stopTime-startTime);

//     // Part c: create points.ppm file in which you draw a bold black circle of radius 3 for every point you generated and you draw a bold red circle of radius 3 for the 2 closest points you found (you may make the circles bold by drawig a circle of radius 2 and a circle of radius 3) 
//     for (auto pt=points.cbegin(); pt!=points.cend(); ++pt) {
//         if ((pt->get_x()==min_pt1.get_x() && pt->get_y()==min_pt1.get_y()) || (pt->get_x()==min_pt2.get_x() && pt->get_y()==min_pt2.get_y())) {
//             Circle c(pt->get_x()*scale_row, pt->get_y()*scale_column, 3);
//             c.draw_circle_color(pixel2D,255,0,0);
//             Circle c2(pt->get_x()*scale_row, pt->get_y()*scale_column, 2);
//             c2.draw_circle_color(pixel2D,255,0,0);
//         } else{
//             Circle c(pt->get_x()*scale_row, pt->get_y()*scale_column, 3);
//             c.draw_circle(pixel2D);
//             Circle c2(pt->get_x()*scale_row, pt->get_y()*scale_column, 2);
//             c2.draw_circle(pixel2D);
//         }
//     }

//     RGB temprgb;
//     for (int r=0; r<scale_row; r++) {
//         for (int c=0; c<scale_column; c++) {
//             temprgb=pixel2D[r*scale_column+c];
//             outfile << temprgb.get_red() << " " << temprgb.get_green() << " " << temprgb.get_blue() << " ";
//         }
//         outfile << "\n";
//     }
//     outfile.close();

//     part1Pair=Pair(min_pt1, min_pt2, min_distance);
//     part1Time=elapsedTime;

//     // Part d: display on the screen the 2 closest points and the distance between them (so the minimum distance)
//     cout << setprecision(17);
//     cout << "The closest pair of points are: \n";
//     cout << "Point 1: (" << min_pt1.get_x() << ", " << min_pt1.get_y() << ")\n";
//     cout << "Point 2: (" << min_pt2.get_x() << ", " << min_pt2.get_y() << ")\n";
//     cout << "The distance between them is: " << min_distance << "\n";
//     cout << "The time for Part One is: " << elapsedTime.count() << " microseconds\n";

//     return 0;
// }

// part 2 code
bool compareX(Point a, Point b) {
    return a.get_x() < b.get_x();
}

bool compareY(Point a, Point b) {
    return a.get_y() < b.get_y();
}

Pair recurClosestPair(vector<Point>& points, int start, int end) {
    // brute force base case
    if (end-start==1) {
        Point A=points[start], B=points[end];
        double minDistance=distance(A, B);
        return Pair(A, B, minDistance);
    } else if (end-start==2) {
        Point A=points[start], B=points[start+1], C=points[end];
        double distAB=distance(A, B), distAC=distance(A, C), distBC=distance(B, C);
        if (distAB==min(distAB, min(distBC,distAC))) { return Pair(A, B, distAB); }
        else if (distAC==min(distAB, min(distBC,distAC))) { return Pair(A, C, distAC); }
        else { return Pair(B, C, distBC); }
    }

    // recursive case
    Pair recurRight=recurClosestPair(points, start, (end+start)/2);
    // cout << "hello there!";
    Pair recurLeft=recurClosestPair(points, (end+start)/2+1, end);
    // cout << "hello there2!";
    Point min_pt1, min_pt2;
    double min;

    // cout << "hello there3!";

    if (recurRight.get_distance()<recurLeft.get_distance()) {
        min_pt1=recurRight.get_pt1();
        min_pt2=recurRight.get_pt2();
        min=recurRight.get_distance();
    } else {
        min_pt1=recurLeft.get_pt1();
        min_pt2=recurLeft.get_pt2();
        min=recurLeft.get_distance();
    }

    Pair result=Pair(min_pt1, min_pt2, min);

    // create the strip
    int midpointIndex=(start+end)/2;
    Point midpoint=points[midpointIndex];

    double leftBound=midpoint.get_x()-min, rightBound=midpoint.get_x()+min;
    int leftBoundIdx=-1, rightBoundIdx=-1;
    bool leftBoundFlag=true, rightBoundFlag=true;

    for (int i=start; i<=midpointIndex; i++) {
        Point temp=points[i];
        if (temp.get_x()>=leftBound && leftBoundFlag==true) {
            leftBoundIdx=i;
            leftBoundFlag=false;
        }
    }

    for (int i=end; i>=midpointIndex+1; i--) {
        Point temp=points[i];
        if (temp.get_x()<=rightBound && rightBoundFlag==true) {
            rightBoundIdx=i;
            rightBoundFlag=false;
        }
    }
    
    if (rightBoundIdx>-1) {
        for (int i=leftBoundIdx; i<=midpointIndex; i++) {
            for (int j=midpointIndex+1; j<=rightBoundIdx; j++) {
                Point pointI=points[i], pointJ=points[j];
                double tempdist=distance(pointI, pointJ);
                if (tempdist<min && tempdist!=0) {
                    min=tempdist;
                    result=Pair(pointI, pointJ, min);
                }
            }
        }
    }
    return result;
}

int part2() {
    vector<Point> points;
    ifstream file("points.txt");
    if (file.fail()) return -1;
    
    // read the points
    double x, y;
    while (file >> x >> y) {
        points.push_back(Point(x, y));
    }
    file.close();
    vector<Point> Px=points, Py=points;
    
    // start timing before sort
    auto startTime2=high_resolution_clock::now();

    // sort points
    sort(Px.begin(), Px.end(), compareX);
    
    Pair finalpair=Pair();
    finalpair=recurClosestPair(Px, 0, points.size()-1);
    // cout << "The closest pair of points are: \n";
    Point min_pt1=finalpair.get_pt1(), min_pt2=finalpair.get_pt2();
    double min_distance=finalpair.get_distance();

    auto stopTime2=high_resolution_clock::now();
    auto elapsedTime=duration_cast<microseconds>(stopTime2-startTime2);

    part2Pair=finalpair;
    part2Time=elapsedTime;

    cout << setprecision(17);
    cout << "The closest pair of points are: \n";
    cout << "Point 1: (" << min_pt1.get_x() << ", " << min_pt1.get_y() << ")\n";
    cout << "Point 2: (" << min_pt2.get_x() << ", " << min_pt2.get_y() << ")\n";
    cout << "The distance between them is: " << min_distance << "\n";
    cout << "The time for Part Two is: " << elapsedTime.count() << " microseconds\n";
    return 0;
}

// part 3 code
Pair optimizedRecurClosestPair(vector<Point>& points, int start, int end) {
    // brute force base case
    if (end-start==1) {
        Point A=points[start], B=points[end];
        double minDistance=distance(A, B);
        return Pair(A, B, minDistance);
    } else if (end-start==2) {
        Point A=points[start], B=points[start+1], C=points[end];
        double distAB=distance(A, B), distAC=distance(A, C), distBC=distance(B, C);
        if (distAB==min(distAB, min(distBC,distAC))) { return Pair(A, B, distAB); }
        else if (distAC==min(distAB, min(distBC,distAC))) { return Pair(A, C, distAC); }
        else { return Pair(B, C, distBC); }
    }
    
    
    // recursive case
    int midpointIndex=(start+end)/2;
    Pair recurRight=optimizedRecurClosestPair(points, start, midpointIndex);
//     cout << "hello there!";
    Pair recurLeft=optimizedRecurClosestPair(points, midpointIndex+1, end);
//     cout << "hello there2!";
    Point min_pt1, min_pt2;
    double minDist;

//     cout << "hello there3!";

    if (recurRight.get_distance()<recurLeft.get_distance()) {
        min_pt1=recurRight.get_pt1();
        min_pt2=recurRight.get_pt2();
        minDist=recurRight.get_distance();
    } else {
        min_pt1=recurLeft.get_pt1();
        min_pt2=recurLeft.get_pt2();
        minDist=recurLeft.get_distance();
    }
    
//     cout << "reached me";
    
    int leftIndex=midpointIndex, rightIndex=midpointIndex+1;
    vector<Point> strip;
    while (points[leftIndex].get_x()>=points[midpointIndex].get_x()-minDist && leftIndex>=start) {
        strip.push_back(points[leftIndex]);
        leftIndex--;
    }
    
    while (points[rightIndex].get_x()<=points[midpointIndex].get_x()+minDist && rightIndex<=end) {
        strip.push_back(points[rightIndex]);
        rightIndex++;
    }
    
    double recurMinDist=minDist;
    sort(strip.begin(), strip.end(), compareY);
    int stripsize=strip.size();
    
    for (int left=0; left<stripsize; left++) {
        for (int right=left+1; right<min(left+16, stripsize); right++) {
            double tempDist=distance(strip[left],strip[right]);
            if ((strip[right].get_y()-strip[left].get_y())>minDist) {
                break;
            }
            if (tempDist>0.0 && tempDist<recurMinDist) {
                recurMinDist=tempDist;
                min_pt1=strip[left];
                min_pt2=strip[right];
            }
        }
    }
//     cout << "I returned " << start << end;
    return Pair(min_pt1, min_pt2, recurMinDist);
}

int part3() {
    vector<Point> points;
    ifstream file("points.txt");
    if (file.fail()) return -1;
    
    // read the points
    double x, y;
    while (file >> x >> y) {
        points.push_back(Point(x, y));
    }
    file.close();
    vector<Point> Px=points, Py=points;
    
    // start timing before sort
    auto startTime3=high_resolution_clock::now();

    // sort points
    sort(Px.begin(), Px.end(), compareX);
    // sort(Py.begin(), Py.end(), compareY);
    
    Pair finalpair=Pair();
    finalpair=optimizedRecurClosestPair(Px, 0, points.size()-1);
    // cout << "The closest pair of points are: \n";
    Point min_pt1=finalpair.get_pt1(), min_pt2=finalpair.get_pt2();
    double min_distance=finalpair.get_distance();

    auto stopTime3=high_resolution_clock::now();
    auto elapsedTime3=duration_cast<microseconds>(stopTime3-startTime3);

    part3Pair=finalpair;
    part3Time=elapsedTime3;

    cout << setprecision(17);
    cout << "The closest pair of points are: \n";
    cout << "Point 1: (" << min_pt1.get_x() << ", " << min_pt1.get_y() << ")\n";
    cout << "Point 2: (" << min_pt2.get_x() << ", " << min_pt2.get_y() << ")\n";
    cout << "The distance between them is: " << min_distance << "\n";
    cout << "The time for Part Three is: " << elapsedTime3.count() << " microseconds\n";
    return 0;
    
//     if you reach a point further than mindist you can just break
}

// part 4 code
class MyHashFunction {
    public:
        size_t operator()(Square s) const { return (hash<unsigned long long>()(s.get_r()))^(hash<unsigned long long>()(s.get_c())); }
};

int part4() {
    vector<Point> points;
    ifstream file("points.txt");
    if (file.fail()) return -1;

    // read the points
    double x, y;
    while (file >> x >> y) {
        points.push_back(Point(x, y));
    }
    file.close();

    // start timing before sort
    auto startTime4=high_resolution_clock::now();
    
    // knuth shuffle points
    int pointssize=points.size();
    srand((unsigned int)time(NULL));
    for (int i=0; i<pointssize; i++) {
        int j=rand()%(pointssize-i)+i ;
        Point temp=points[i];
        points[i]=points[j];
        points[j]=temp;
    }

    // map setup
    unordered_map<Square, Point, MyHashFunction> map;
    double min_distance=distance(points[0], points[1]);
    double s=min_distance/2.0;
    Pair finalpair=Pair(points[0], points[1], min_distance);

    // search squares
    for (unsigned long long i=0; i<points.size(); i++) {
        Point temp=points[i];
        unsigned long long temp_row=(unsigned long long)(floor(temp.get_x()/s));
        unsigned long long temp_col=(unsigned long long)(floor(temp.get_y()/s));
        Square temp_square=Square(temp_row, temp_col);

        if (i<=1) { map[temp_square]=temp; }

        unsigned long long rowbounds=temp_row<2?0:temp_row-2;
        unsigned long long colbounds=temp_col<2?0:temp_col-2;

        bool flag=false;
        for (unsigned long long r=rowbounds; r<=temp_row+2; r++) {
            for (unsigned long long c=colbounds; c<=temp_col+2; c++) {
                Square key=Square(r, c);
                if (map.find(key)!=map.end()) {
                    Point squarepoint=map.at(key);
                    double temp_distance=distance(temp, squarepoint);
                    if (temp_distance<min_distance&&temp_distance>0) {
                        min_distance=temp_distance;
                        finalpair=Pair(temp, squarepoint, min_distance);
                        flag=true;
                    }
                }
            }
        }

        if (flag) {
            s=min_distance/2.0;
            map.clear();
            for (unsigned long long x=0; x<=i; x++) {
                Point pasttemp=points[x];
                unsigned long long pasttemp_row=(unsigned long long)(floor(pasttemp.get_x() / s));
                unsigned long long pasttemp_col=(unsigned long long)(floor(pasttemp.get_y() / s));
                Square past_square=Square(pasttemp_row, pasttemp_col);
                map[past_square]=pasttemp;
            }
        } else {
            map[temp_square]=temp;
        }
    }
    
    Point min_pt1=finalpair.get_pt1(), min_pt2=finalpair.get_pt2();
    
    auto stopTime4=high_resolution_clock::now();
    auto elapsedTime4=duration_cast<microseconds>(stopTime4-startTime4);

    part4Pair=finalpair;
    part4Time=elapsedTime4;

    cout << setprecision(17);
    cout << "The closest pair of points are: \n";
    cout << "Point 1: (" << min_pt1.get_x() << ", " << min_pt1.get_y() << ")\n";
    cout << "Point 2: (" << min_pt2.get_x() << ", " << min_pt2.get_y() << ")\n";
    cout << "The distance between them is: " << finalpair.get_distance() << "\n";
    cout << "The time for Part Four is: " << elapsedTime4.count() << " microseconds\n";
    return 0;
}

int main() {
    part0();
    part3();
    part4();
    
    ofstream outfile;
    outfile.open("results.txt");
    outfile << setprecision(17);
    // outfile << "The closest pair of points for PART 1 is: \n";
    // outfile << "Point 1: (" << part1Pair.get_pt1().get_x() << ", " << part1Pair.get_pt1().get_y() << ")\n";
    // outfile << "Point 2: (" << part1Pair.get_pt2().get_x() << ", " << part1Pair.get_pt2().get_y() << ")\n";
    // outfile << "The distance between them is: " << part1Pair.get_distance() << "\n";
    // outfile << "The time is: " << part1Time.count() << " microseconds\n\n"; 

    // outfile << "The closest pair of points for PART 2 is: \n";
    // outfile << "Point 1: (" << part2Pair.get_pt1().get_x() << ", " << part2Pair.get_pt1().get_y() << ")\n";
    // outfile << "Point 2: (" << part2Pair.get_pt2().get_x() << ", " << part2Pair.get_pt2().get_y() << ")\n";
    // outfile << "The distance between them is: " << part2Pair.get_distance() << "\n";
    // outfile << "The time is: " << part2Time.count() << " microseconds\n\n"; 

    outfile << "The closest pair of points for PART 3 is: \n";
    outfile << "Point 1: (" << part3Pair.get_pt1().get_x() << ", " << part3Pair.get_pt1().get_y() << ")\n";
    outfile << "Point 2: (" << part3Pair.get_pt2().get_x() << ", " << part3Pair.get_pt2().get_y() << ")\n";
    outfile << "The distance between them is: " << part3Pair.get_distance() << "\n";
    outfile << "The time is: " << part3Time.count() << " microseconds\n\n"; 

    outfile << "The closest pair of points for PART 4 is: \n";
    outfile << "Point 1: (" << part4Pair.get_pt1().get_x() << ", " << part4Pair.get_pt1().get_y() << ")\n";
    outfile << "Point 2: (" << part4Pair.get_pt2().get_x() << ", " << part4Pair.get_pt2().get_y() << ")\n";
    outfile << "The distance between them is: " << part4Pair.get_distance() << "\n";
    outfile << "The time is: " << part4Time.count() << " microseconds\n\n"; 
    outfile.close();
}