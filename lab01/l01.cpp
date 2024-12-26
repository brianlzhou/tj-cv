#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <limits>
using namespace std;

// Pixel-related commands

// auto pixel2D = new bool[800][800];
auto pixel2D = new bool[800][800];
void fill_pixel(int x, int y) {
    if(x<0 || x>=800 || y<0 || y>=800) return;
    pixel2D[x][y] = true;
}

// Basic functions

double calc_slope(pair<double, double> a, pair<double, double> b) {
    double xA = a.first, yA = a.second;
    double xB = b.first, yB = b.second;

    if (xB==xA) return std::numeric_limits<double>::infinity();
    return (yB-yA)/(xB-xA);
}

double calc_distance(pair<double, double> a, pair<double, double> b) {
    double xA = a.first, yA = a.second;
    double xB = b.first, yB = b.second;

    return sqrt(pow(yB-yA, 2)+pow(xB - xA, 2));
}

pair<double, double> calc_midpoint(pair<double, double> a, pair<double, double> b) {
    double xA = a.first, yA = a.second;
    double xB = b.first, yB = b.second;

    return make_pair((xA+xB)/2, (yA+yB)/2);
}

pair<double, double> rand_point() {
    return make_pair(((double) rand()/(RAND_MAX+1.0)), ((double) rand()/(RAND_MAX+1.0)));
}

// Bressenham's Algorithm

void draw_line_xP(pair<int, int> a, pair<int, int> b){
    int dX=b.first-a.first,
        dY=b.second-a.second,
        y=a.second,
        dYdXdiff=dY-dX;
    for (int x=a.first; x<=b.first-1; x++) {
        fill_pixel(x, y);
        if (dYdXdiff>0) {
            y += 1;
            dYdXdiff -= dX;
        }
        dYdXdiff+=dY;
    }
    int mx = 0;
}

void draw_line_xN(pair<int, int> a, pair<int, int> b) {
    int dX=b.first-a.first,
        dY=b.second-a.second,
        y=a.second,
        dYdXdiff=dX-dY;
    for (int x=a.first; x<=b.first-1; x++) {
        fill_pixel(x, y);
        if (dYdXdiff>0) {
            y-=1;
            dYdXdiff-=dX;
        }
        dYdXdiff+=abs(dY);
    }
}

void draw_line_yP(pair<int, int> a, pair<int, int> b) {
    int dX = b.first - a.first, 
        dY = b.second - a.second, 
        x = a.first, 
        dYdXdiff = dX - dY;
    for (int y=a.second; y<=b.second-1; y++) {
        fill_pixel(x, y);
        if (dYdXdiff > 0) {
            x+=1;
            dYdXdiff-=dY;
        }
        dYdXdiff+=dX;
    }
}

void draw_line_yN(pair<int, int> a, pair<int, int> b) {
    int dX = b.first - a.first, 
        dY = b.second - a.second, 
        x = a.first, 
        dYdXdiff = dX - dY;
    for (int y=a.second; y<=b.second-1; y++) {
        fill_pixel(x, y);
        if (dYdXdiff>0) {
            x-=1;
            dYdXdiff-=dY;
        }
        dYdXdiff+=abs(dX);
    }
}

void draw_line(pair<double, double> a, pair<double, double> b) {
    pair<int, int> scaledA = make_pair((int)(a.first * 800), (int)(a.second * 800));
    pair<int, int> scaledB = make_pair((int)(b.first * 800), (int)(b.second * 800));

    int dx=scaledB.first-scaledA.first;
    int dy=scaledB.second-scaledA.second;

    bool isAbsDxGreater = abs(dx) > abs(dy);
    
    if(dx>=dy) 
        if(dy<0 && dx<0) draw_line_yP(scaledB, scaledA);
        else if(dy < 0) 
            if(abs(dx) > abs(dy))
                draw_line_xN(scaledA, scaledB);
            else
                draw_line_yN(scaledB, scaledA);
        else draw_line_xP(scaledA, scaledB);
    else if(dy<0 && dx<0) draw_line_xP(scaledB, scaledA);
    else if(dx < 0) 
        if(abs(dx) < abs(dy))
            draw_line_yN(scaledA, scaledB);
        else
            draw_line_xN(scaledB, scaledA);
    else draw_line_yP(scaledA, scaledB);
}

// Circle functions

void draw_circle(pair<double, double> center, double r) {
    pair<int, int> scaledCenter = make_pair((int)(center.first * 800), (int)(center.second * 800));
    int scaledR = (int)(r * 800);

    int x, y, xmax, y2, y2_new, ty;

    xmax = (int) (scaledR * 0.70710678);

    y=scaledR;
    y2 = y * y;
    ty = (2 * y) - 1;
    y2_new = y2;

    for (int x = 0; x <= xmax+1; x++) {
        if ((y2 - y2_new) >= ty) {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        fill_pixel(scaledCenter.first+x, scaledCenter.second-y);
        fill_pixel(scaledCenter.first-x, scaledCenter.second-y);
        fill_pixel(scaledCenter.first+x, scaledCenter.second+y);
        fill_pixel(scaledCenter.first-x, scaledCenter.second+y);
        fill_pixel(scaledCenter.first+y, scaledCenter.second-x);
        fill_pixel(scaledCenter.first-y, scaledCenter.second-x);
        fill_pixel(scaledCenter.first+y, scaledCenter.second+x);
        fill_pixel(scaledCenter.first-y, scaledCenter.second+x);
        y2_new -= (2 * x) - 3;
    }   
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an alternate Bressenham's circle algorithm implementation not used in the code or main function.
// //////////////////////////////////////////////////////////////////////////////////////////////////////////
void draw_circle_better(pair<double, double> center, double r) {
    pair<int, int> scaledCenter = make_pair((int)(center.first * 800), (int)(center.second * 800));
    int scaledR = (int)(r * 800);
    int x = scaledR, y = 0;
    
    // Function to plot eight symmetric points for the given point
    auto plot_symmetric_points = [&](int x, int y) {
        fill_pixel(scaledCenter.first+x, scaledCenter.second-y);
        fill_pixel(scaledCenter.first-x, scaledCenter.second-y);
        fill_pixel(scaledCenter.first+x, scaledCenter.second+y);
        fill_pixel(scaledCenter.first-x, scaledCenter.second+y);
        fill_pixel(scaledCenter.first+y, scaledCenter.second-x);
        fill_pixel(scaledCenter.first-y, scaledCenter.second-x);
        fill_pixel(scaledCenter.first+y, scaledCenter.second+x);
        fill_pixel(scaledCenter.first-y, scaledCenter.second+x);
    };

    plot_symmetric_points(x, y);

    int decisionParam = 1 - scaledR;
    while (x>y) {
        y++;
        if 
            (decisionParam<=0) decisionParam=decisionParam+2*y+1;
        else{
            x--;
            decisionParam=decisionParam+2*y-2*x+1;
        } 
        plot_symmetric_points(x, y);
    }
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an alternate Bressenham's circle algorithm implementation not used in the code or main function.
// //////////////////////////////////////////////////////////////////////////////////////////////////////////

pair<double, double> calc_circumradius_inradius(double a, double b, double c)
{
    double s = (a+b+c)/2; // semiperimeter

    double circumradius = sqrt((s-a)*(s-b)*(s-c)/s);
    double inradius = (a*b*c)/(4*circumradius*s);

    return make_pair(circumradius, inradius);
}

void draw_ninepoint_circle(pair<double, double> point1, 
                           pair<double, double> point2, 
                           pair<double, double> point3,
                           double centerIncircleX, double centerIncircleY, 
                           const pair<double, double> bisectorC,
                           const pair<double, double> bisectorB) {
    double originalSlope1=calc_slope(point2, point1);
    double originalSlope2=calc_slope(point3, point1);
    double originalSlope3=calc_slope(point3, point2);
    double perpendicularSlope1=-1/originalSlope1;
    double perpendicularSlope2=-1/originalSlope2;
    double perpendicularSlope3=-1/originalSlope3;

    double orthogonalX = (point2.second-perpendicularSlope2*point2.first+perpendicularSlope3*point1.first-point1.second)/(perpendicularSlope3 - perpendicularSlope2);
    double orthogonalY = point1.second + perpendicularSlope3 * (orthogonalX - point1.first);

    // Draw the euler line
    double slopeCircOrtho = calc_slope(make_pair(centerIncircleX, centerIncircleY), make_pair(orthogonalX, orthogonalY));
    double intercept = centerIncircleY - (slopeCircOrtho * centerIncircleX);
    draw_line(make_pair(0, intercept), make_pair(1, slopeCircOrtho+intercept));

    // Draw the nine-point circle
    auto centerNinecircle=calc_midpoint(make_pair(centerIncircleX, centerIncircleY), make_pair(orthogonalX, orthogonalY));
    double radiusNinecircle=calc_distance(centerNinecircle, calc_midpoint(point1, point2));
    draw_circle(centerNinecircle, radiusNinecircle);
    cout << "Nine-point circle | radius: " << radiusNinecircle << " | center: (" << centerNinecircle.first << ", " << centerNinecircle.second << ")" << endl;
    cout << "Euler line | slope: " << slopeCircOrtho << " | y-intercept: " << intercept << endl;
}

pair<double, double> calc_perpendicular_bisector(pair<double, double> pt1, pair<double, double> pt2) {
    double originalSlope = calc_slope(pt1, pt2);
    double perpendicularSlope = -1.0/originalSlope;

    auto midpoint = calc_midpoint(pt1, pt2);
    double yIntercept = midpoint.second-(perpendicularSlope*midpoint.first);

    return make_pair(perpendicularSlope, yIntercept);
}

int main() {
    ofstream outfile;
    outfile.open("triangles.ppm");
    srand(time(0));
    outfile << "P3\n" << 800 << " " << "800\n" << "1\n";

    auto point1=rand_point();
    auto point2=rand_point();
    auto point3=rand_point();
    // // test case 1
    // point1=make_pair(0.563568, 0.00125122) ;
    // point2=make_pair(0.808716, 0.193298) ;
    // point3=make_pair(0.479858, 0.584991) ;
    
    // // test case 2
    // point1=make_pair(.24,.75);
    // point2=make_pair(.1,.5);
    // point3=make_pair(.7,.23);
    cout << "Point 1: (" << point1.first << ", " << point1.second << ") (" << 800*point1.first << ", " << 800*point1.second << ")" << endl;
    cout << "Point 2: (" << point2.first << ", " << point2.second << ") (" << 800*point2.first << ", " << 800*point2.second << ")" << endl;
    cout << "Point 3: (" << point3.first << ", " << point3.second << ") (" << 800*point3.first << ", " << 800*point3.second << ")" << endl;

    double pt1X=point1.first, pt1Y=point1.second;
    double pt2X=point2.first, pt2Y=point2.second;
    double pt3X=point3.first, pt3Y=point3.second;

    double edgeA=calc_distance(point1, point2);
    double edgeB=calc_distance(point2, point3);
    double edgeC=calc_distance(point1, point3);

    // Draw the triangle
    draw_line(point1, point2);
    draw_line(point2, point3); 
    draw_line(point1, point3);

    // Draw the test circle
    draw_circle(make_pair(.25,.75), 0.00625);

    // Calculate incircle, circumcircle radius
    auto ciurcuminradii = calc_circumradius_inradius(edgeA, edgeB, edgeC);
    
    // Draw the incircle
    auto centerIncircle = make_pair((pt1X*edgeB+pt2X*edgeC+pt3X*edgeA)/(edgeA+edgeB+edgeC), (pt1Y*edgeB+pt2Y*edgeC+pt3Y*edgeA)/(edgeA+edgeB+edgeC));
    draw_circle(centerIncircle, ciurcuminradii.first);
    cout << "Incircle | radius: " << ciurcuminradii.first << " | center: (" << centerIncircle.first << ", " << centerIncircle.second << ")" << endl;

    // Draw the circumcircle
    auto bisectorC = calc_perpendicular_bisector(point1, point2);
    auto bisectorA = calc_perpendicular_bisector(point2, point3);
    auto bisectorB = calc_perpendicular_bisector(point1, point3);
    double centerIncircleX=(bisectorA.second-bisectorC.second)/(bisectorC.first-bisectorA.first);
    double centerIncircleY=(centerIncircleX*bisectorC.first) + bisectorC.second;
    draw_circle(make_pair(centerIncircleX,centerIncircleY), ciurcuminradii.second);
    cout << "Circumcircle | radius: " << ciurcuminradii.second << " | center: (" << centerIncircleX << ", " << centerIncircleY << ")" << endl;

    // Draw the 9-point circle & Euler line
    draw_ninepoint_circle(point1, point2, point3, centerIncircleX, centerIncircleY, bisectorC, bisectorB);

    // Write ppm
    for(int i=0; i<800; i++) {
        for(int j=0; j<800; j++) {
            if(pixel2D[i][j]==true) outfile << 0 << " " << 0 << " " << 0 << " ";
            else outfile << 1 << " " << 1 << " " << 1 << " ";
        }
    }
    outfile.close();
    std::cout << "successfully generated triangles.ppm!" << std::endl;

    return 0;
}