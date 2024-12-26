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

int window_row=0;
int window_column=0;

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

RGB* pixel2D = new RGB[window_row*window_column];

void fill_pixel(RGB* imgcanvas, int x, int y) {
    if(x<0 || x>=400 || y<0 || y>=400) return;
    *(imgcanvas+y*window_row+x)=RGB(0,0,0);
}

void fill_pixel_color(RGB* imgcanvas, int x, int y, int r, int g, int b) {
    if(x<0 || x>=400 || y<0 || y>=400) return;
    *(imgcanvas+y*window_row+x)=RGB(r,g,b);
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
                        *(imgcanvas+idx*window_row+x1) = RGB(0,0,0);
                    }
                } else {
                    for (int idx=y2; idx<=y1; idx++) {
                        *(imgcanvas+idx*window_row+x1) = RGB(0,0,0);
                    }
                }
            } else {
                if (x2>x1) {
                    for (int idx=x1; idx<=x2; idx++) {
                        *(imgcanvas+y1*window_row+idx) = RGB(0,0,0);
                    }
                } else {
                    for (int idx=x2; idx<=x1; idx++) {
                        *(imgcanvas+y1*window_row+idx) = RGB(0,0,0);
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
                *(imgcanvas+y*window_row+x)=RGB(0,0,0);
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
                *(imgcanvas+y*window_row+x)=RGB(0,0,0);
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
                *(imgcanvas+y*window_row+x)=RGB(0,0,0);
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
                *(imgcanvas+y*window_row+x)=RGB(0,0,0);
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

int part1() {
    string p3or6;
    int rgb;

    ifstream imagefile("image.ppm");
    imagefile >> p3or6 >> window_column >> window_row >> rgb;

    // grayscale
    vector<vector<int>> grayscale(window_row, vector<int>(window_column));
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            int red,green,blue;
            imagefile >> red >> green >> blue;
            int weightedAverage=round(0.299*red+0.587*green+0.114*blue); 
            grayscale[r][c]=weightedAverage;
        }
    }
    imagefile.close();
    
    // sobel
    vector<vector<int>> sobel(window_row, vector<int>(window_column));
    int single_threshold=150;
    for (int r=1; r<window_row-1; r++) {
        for (int c=1; c<window_column-1; c++) {
            int s_grad=round(sqrt(pow(grayscale[r-1][c+1]-grayscale[r+1][c-1] + 2*grayscale[r][c+1]-2*grayscale[r][c-1] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)+
                                  pow(grayscale[r+1][c-1]-grayscale[r-1][c+1] + 2*grayscale[r+1][c]-2*grayscale[r-1][c] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)));
            if (s_grad>single_threshold) { sobel[r][c]=1; }
            else{ sobel[r][c]=0; }
        }
    }
    
    // write
    ofstream outfile;
    outfile.open("imageg.ppm");
    outfile << "P3 " << window_column << " " << window_row << " 255\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfile << grayscale[r][c] << " " << grayscale[r][c] << " " << grayscale[r][c] << " ";
        }
        outfile << "\n";
    }
    outfile.close();
    
    ofstream outfile_sobel;
    outfile_sobel.open("imagem.ppm");
    outfile_sobel << "P3 " << window_column << " " << window_row << " 1\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfile_sobel << sobel[r][c] << " " << sobel[r][c] << " " << sobel[r][c] << " ";
        }
        outfile_sobel << "\n";
    }
    outfile_sobel.close();
    return 0;
}

int main() {
    part1();
}