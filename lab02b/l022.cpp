#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <algorithm>
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
};

// Part 1 code
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

// Part 2 code
Point find_intersection(Point pt1, Point pt2, double slope1, double slope2, bool is_vertical) {
    // break down point 1 and 2 into x and y components
    double pt1x=pt1.get_x(), pt1y=pt1.get_y();
    double pt2x=pt2.get_x(), pt2y=pt2.get_y();

    // calculate the intersection point
    double intersectionX, intersectionY;
    if (is_vertical) { // if the line is vertical (must go first; slope is undef)
        intersectionX=pt2x;
        intersectionY=pt1y;
        return Point(intersectionX, intersectionY);
    }
    else if (slope2==0) {
        intersectionX=pt1x;
        intersectionY=pt2y;
        return Point(intersectionX, intersectionY);
    }
    else {
        double a1, b1, a2, b2, c1, c2, dt, dx, dy;
        a1=-slope1, a2=-slope2;
        b1=1, b2=1;
        c1=pt1y-slope1*pt1x, c2=pt2y-slope2*pt2x;

        dt=a1*b2-a2*b1;
        dx=c1*b2-c2*b1;
        dy=a1*c2-a2*c1;

        intersectionX=dx/dt;
        intersectionY=dy/dt;
        return Point(intersectionX, intersectionY);
    }
}

void calculate_various_slopes(Point square1, Point square3, Point square2, Point square4, Point equidistantpt, 
                                     Point squares[square_combinations][square_pts], double areas[], int n) {
//     double pt1x=square1.get_x(), pt1y=square1.get_y();
    double pt2x=square2.get_x(), pt2y=square2.get_y();
//     double pt3x=square3.get_x(), pt3y=square3.get_y();
//     double pt4x=square4.get_x(), pt4y=square4.get_y();
    double ptEx=equidistantpt.get_x(), ptEy=equidistantpt.get_y();

    // Step 5: the last point that you didn't use and E are one of the lines of the square
    // calculate slopes of various lines
    double slopefor_be; 
    bool vertical_be=false;
    if (ptEx==pt2x) vertical_be=true;
    else slopefor_be=(pt2y-ptEy)/(pt2x-ptEx);

    double slopefor_ac; 
    if (vertical_be) slopefor_ac=0;
    else if (slopefor_be!=0) slopefor_ac=-1.0/slopefor_be;

    double slopefor_d; 
    bool vertical_d=false;
    if (vertical_be) vertical_d= true;
    else slopefor_d=slopefor_be;

    // Step 6: from the 2 points that you pick at step 1, draw perpendiculars on the line from step 5
    // calculate the intersection point of various lines
    Point intersectionAB, intersectionBC, intersectionCD, intersectionDA;
    intersectionAB = find_intersection(square1, square2, slopefor_ac, slopefor_be, vertical_be);
    intersectionBC = find_intersection(square3, square2, slopefor_ac, slopefor_be, vertical_be);
    intersectionCD = find_intersection(square3, square4, slopefor_ac, slopefor_d, vertical_d);
    intersectionDA = find_intersection(square1, square4, slopefor_ac, slopefor_d, vertical_d);
    
    double area=pow(intersectionAB.get_x()-intersectionBC.get_x(), 2)+pow(intersectionAB.get_y()-intersectionBC.get_y(), 2);

    squares[n][0]=intersectionAB;
    squares[n][1]=intersectionBC;
    squares[n][2]=intersectionCD;
    squares[n][3]=intersectionDA;
    areas[n] = area;
}

void calc_square(Point square1, Point square3, Point square2, Point square4, Point squares[square_combinations][square_pts], double areas[], int n) {
    double pt1x=square1.get_x(), pt1y=square1.get_y();
//     double pt2x=square2.get_x(), pt2y=square2.get_y();
    double pt3x=square3.get_x(), pt3y=square3.get_y();
    double pt4x=square4.get_x(), pt4y=square4.get_y();

    // Step 3: draw a perpendicular line from the point in step 2 to the line in step 1
    // calculate the slope of the line and perpendicular line
    double opposite_point_slope;
    bool opposite_points_are_vertical = false;
    double opposite_point_perp_slope;
    bool opposite_points_are_horizontal = false; //same as if perpendicular line is vertical

    if (pt1x==pt3x) opposite_points_are_vertical = true;
    else opposite_point_slope = (pt3y-pt1y)/(pt3x-pt1x);

    if (opposite_points_are_vertical) opposite_point_perp_slope=0;
    else if (opposite_point_slope==0) opposite_points_are_horizontal=true;
    else opposite_point_perp_slope=-1.0/opposite_point_slope;

    // calculate the euclidian distnace between the opposite points
    double euclidian_btwn13 = sqrt(pow((pt1x-pt3x), 2)+pow((pt1y-pt3y), 2));
    
    // Step 4: pick E a point on the line from step 4 such that the point from step 3 and E is equal to the segment from step 1
    // calculate the points that go in either direction
    Point equidistant1, equidistant2;
    if (opposite_points_are_horizontal) {
        equidistant1=Point(pt4x, pt4y+euclidian_btwn13);
        equidistant2=Point(pt4x, pt4y-euclidian_btwn13);
    } else {
        double tempx=pt4x+sqrt(pow(euclidian_btwn13,2)/(pow(opposite_point_perp_slope,2)+1));
        equidistant1=Point(tempx, pt4y+opposite_point_perp_slope*(tempx-pt4x));
        double tempx2=pt4x-sqrt(pow(euclidian_btwn13,2)/(pow(opposite_point_perp_slope,2)+1));
        equidistant2=Point(tempx2, pt4y+opposite_point_perp_slope*(tempx2-pt4x));
    }

    // calculate both of the squares
    calculate_various_slopes(square1, square3, square2, square4, equidistant1, squares, areas, (n-1)*2);
    cout <<  setprecision(17);
    // cout << square1.get_x() << square1.get_y() << square3.get_x() << square3.get_y() << square2.get_x() << square2.get_y() << square4.get_x() << square4.get_y() << equidistant1.get_x() << equidistant1.get_y() << endl;
    calculate_various_slopes(square1, square3, square2, square4, equidistant2, squares, areas, (n-1)*2+1);
}

int find_smallest_square(Point squares[square_combinations][square_pts], Point ptA, Point ptB, Point ptC, Point ptD) {   
    // Step 1: pick 2 points and they will be on opposite sides of the square
    // assigned by arranging points in calc_square, first two are the picked points
    // Step 2: pick one of the other points (no effect)
    double areas[6];
    calc_square(ptA, ptB, ptC, ptD, squares, areas, 1);
    calc_square(ptA, ptC, ptB, ptD, squares, areas, 2);
    calc_square(ptA, ptD, ptB, ptC, squares, areas, 3);

    ofstream outputfile;
    outputfile.open("output.txt");
    outputfile << setprecision(17);
    outputfile << "(" << ptA.get_x() << "," << ptA.get_y() << ") , ";
    outputfile << "(" << ptB.get_x() << "," << ptB.get_y() << ") , ";
    outputfile << "(" << ptC.get_x() << "," << ptC.get_y() << ") , ";
    outputfile << "(" << ptD.get_x() << "," << ptD.get_y() << ")\n";

    double smallestarea=10000000000.1;
    int smallestarea_index=0;
    for (int row=0; row<square_combinations; row++) {
        for (int col=0; col<square_pts; col++) {
            Point temppt = squares[row][col];
            if (col<=2) outputfile << "(" << temppt.get_x() << "," << temppt.get_y() << ") , ";
            else outputfile << "(" << temppt.get_x() << "," << temppt.get_y() << ") ";
        }
        outputfile << "Area=" << areas[row] << endl;
        if (areas[row]<smallestarea) {
            smallestarea=areas[row];
            smallestarea_index=row;
        }
    }
    outputfile.close();
    return smallestarea_index;
}

void calc_endpoints(double x1, double y1, double x2, double y2, int* x3, int* y3, int* x4, int* y4) {
    if (x1==x2) {
        *x3=std::round(x1);
        *y3=0;
        *x4=std::round(x1);
        *y4=scale_row-1;
    } 
    else if (y1==y2) {
        *x3=0;
        *y3=round(y1);
        *x4=scale_column-1;
        *y4=round(y1);
    } 
    else {
        double m=(y2-y1)/(x2-x1);
        double b=y1-m*x1;

        int x_at_ymin = round((0-b)/m);
        int x_at_ymax = round((scale_row-1-b)/m);
        int y_at_xmin = round(m*0+b);
        int y_at_xmax = round(m*(scale_column-1)+b);
        
        if (0 <= x_at_ymin && x_at_ymin<scale_column) {
            *x3=x_at_ymin;
            *y3=0;
        } else if (0 <= y_at_xmin && y_at_xmin<scale_row) {
            *x3=0;
            *y3=y_at_xmin;
        }

        if (0<=x_at_ymax && x_at_ymax<scale_column) {
            *x4=x_at_ymax;
            *y4=scale_row-1;
        } else if (0 <= y_at_xmax && y_at_xmax<scale_row) {
            *x4=scale_column-1;
            *y4=y_at_xmax;
        }
    }
}

int part2() {
    // setup output ppm file
    ofstream outfile;
    outfile.open("output.ppm");
    outfile << "P3 " << scale_column << " " << scale_row << " " << "255" << endl;

    RGB* pixel2D=new RGB[scale_column*scale_row];
    for (int c=0; c<scale_column; c++) {
        for (int r=0; r<scale_row; r++) {
            *(pixel2D+r*scale_column+c)=RGB();
        }
    }

    // read in points
    string templine; 
    std::ifstream reader("points.txt");
    if (reader.fail()) return -1;
    getline(reader, templine);
    std::string::iterator iter;
    string new_templine;

    for (iter=templine.begin(); iter!=templine.end(); iter++) {
        if (*iter!=')' && *iter!='(' && *iter!=' ') new_templine.push_back(*iter);
    }

    string points2D[8];
    stringstream ss(new_templine);
    string word;
    int i = 0;
    while (!ss.eof()) {
        getline(ss, word, ',');
        points2D[i]=word;
        i++;
    }

    Point ptA = Point(stod(points2D[0]),stod(points2D[1]));
    Point ptB = Point(stod(points2D[2]),stod(points2D[3]));
    Point ptC = Point(stod(points2D[4]),stod(points2D[5]));
    Point ptD = Point(stod(points2D[6]),stod(points2D[7]));

    // scale all four original points by the 800 scale
    int xA=int(round(scale_column*ptA.get_x())), yA=int(round(scale_row*ptA.get_y()));
    int xB=int(round(scale_column*ptB.get_x())), yB=int(round(scale_row*ptB.get_y()));
    int xC=int(round(scale_column*ptC.get_x())), yC=int(round(scale_row*ptC.get_y()));
    int xD=int(round(scale_column*ptD.get_x())), yD=int(round(scale_row*ptD.get_y()));

    Circle circA=Circle(xA, yA, 2);
    circA.draw_circle(pixel2D);
    Circle circB=Circle(xB, yB, 2);
    circB.draw_circle(pixel2D);
    Circle circC=Circle(xC, yC, 2);
    circC.draw_circle(pixel2D);
    Circle circD=Circle(xD, yD, 2);
    circD.draw_circle(pixel2D);

    Point all_squares[6][4];
    Point temp;
    for (int row=0; row < 6; row++) {
        for (int col = 0; col < 4; col++) {
            all_squares[row][col]=Point(0.0,0.0);
        }
    }

    int smallest_square;
    smallest_square=find_smallest_square(all_squares, ptA, ptB, ptC, ptD);
    Point squareA = all_squares[smallest_square][0];
    Point squareB = all_squares[smallest_square][1];
    Point squareC = all_squares[smallest_square][2];
    Point squareD = all_squares[smallest_square][3];

    // scale all four points of the smallest square by the 800 scale
    int square_xA=int(round(scale_column*squareA.get_x())), square_yA=int(round(scale_row*squareA.get_y()));
    int square_xB=int(round(scale_column*squareB.get_x())), square_yB=int(round(scale_row*squareB.get_y()));
    int square_xC=int(round(scale_column*squareC.get_x())), square_yC=int(round(scale_row*squareC.get_y()));
    int square_xD=int(round(scale_column*squareD.get_x())), square_yD=int(round(scale_row*squareD.get_y()));

    // draw circles around all four points of the smallest square
    Circle circ_squareA=Circle(square_xA, square_yA, 4);
    circ_squareA.draw_circle(pixel2D);
    Circle circ_squareB=Circle(square_xB, square_yB, 4);
    circ_squareB.draw_circle(pixel2D);
    Circle circ_squareC=Circle(square_xC, square_yC, 4);
    circ_squareC.draw_circle(pixel2D);
    Circle circ_squareD=Circle(square_xD, square_yD, 4);
    circ_squareD.draw_circle(pixel2D);

    // draw the smallest square
    double square_DxA=scale_column*squareA.get_x(), square_DyA=scale_row*squareA.get_y();
    double square_DxB=scale_column*squareB.get_x(), square_DyB=scale_row*squareB.get_y();
    double square_DxC=scale_column*squareC.get_x(), square_DyC=scale_row*squareC.get_y();
    double square_DxD=scale_column*squareD.get_x(), square_DyD=scale_row*squareD.get_y();

    int xEND1, yEND1, xEND2, yEND2;
    calc_endpoints(square_DxA, square_DyA, square_DxB, square_DyB, &xEND1, &yEND1, &xEND2, &yEND2);
    LineSegment sideAB = LineSegment(xEND1, yEND1, xEND2, yEND2);
    sideAB.draw_line(pixel2D);
    calc_endpoints(square_DxB, square_DyB, square_DxC, square_DyC, &xEND1, &yEND1, &xEND2, &yEND2);
    LineSegment sideBC = LineSegment(xEND1, yEND1, xEND2, yEND2);
    sideBC.draw_line(pixel2D);
    calc_endpoints(square_DxC, square_DyC, square_DxD, square_DyD, &xEND1, &yEND1, &xEND2, &yEND2);
    LineSegment sideCD = LineSegment(xEND1, yEND1, xEND2, yEND2);
    sideCD.draw_line(pixel2D);
    calc_endpoints(square_DxD, square_DyD, square_DxA, square_DyA, &xEND1, &yEND1, &xEND2, &yEND2);
    LineSegment sideDA = LineSegment(xEND1, yEND1, xEND2, yEND2);
    sideDA.draw_line(pixel2D);

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
    part2();
}