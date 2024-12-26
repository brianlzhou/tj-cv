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
    int single_threshold=80; // value shown in horse.jpg
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

void applyHysteresis(vector<vector<int>>& edges, vector<vector<int>>& visited, int x, int y, int lowThreshold, int highThreshold) {
    // if hit bound, already a strong edge, or less than low threshold:
    if (x<0 || x>=window_row || y<0 || y>=window_column || visited[x][y]==1 || edges[x][y]<lowThreshold) { return; }

    // mark as strong edge:
    edges[x][y]=255;
    visited[x][y]=1;

    // recur to neighbors:
    for (int dx=-1; dx<=1; dx++) {
        for (int dy=-1; dy<=1; dy++) {
            if (dx!=0 || dy!=0) { applyHysteresis(edges, visited, x+dx, y+dy, lowThreshold, highThreshold); }
        }
    }
}

int part2(int argc, char *argv[]) {
    string p3or6;
    int rgb;
    int lowThreshold=240, highThreshold=300; // values 240/300 optimised for horse.jpg
    string infilename=("image.ppm"), outfilename=("image1.ppm");
    
    if (argc>1) {
        for (int i=0; i+2<argc; i+=2) {
            string curr=argv[i+1];
            if (curr=="-f") {infilename=argv[i+2];}
            if (curr=="-lt") {
                string lowThreshStr=argv[i+2];
                lowThreshold=stoi(lowThreshStr);
            }
            if (curr=="-lt") {
                string highThreshStr=argv[i+2];
                highThreshold=stoi(highThreshStr);
            }
            if (curr=="-of") {outfilename=argv[i+2];}
        }
    }
    
    ifstream imagefile(infilename);
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
    for (int r=1; r<window_row-1; r++) {
        for (int c=1; c<window_column-1; c++) {
            int s_grad=round(sqrt(pow(grayscale[r-1][c+1]-grayscale[r+1][c-1] + 2*grayscale[r][c+1]-2*grayscale[r][c-1] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)+
                                  pow(grayscale[r+1][c-1]-grayscale[r-1][c+1] + 2*grayscale[r+1][c]-2*grayscale[r-1][c] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)));
            sobel[r][c] = s_grad;
        }
    }

    // double thresh
    vector<vector<int>> visited(window_row, vector<int>(window_column, 0));
    vector<vector<int>> edges(window_row, vector<int>(window_column, 0));
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (sobel[r][c]>highThreshold) {
                edges[r][c]=255; // Strong edge
            } else if (sobel[r][c]>lowThreshold) {
                edges[r][c]=sobel[r][c]; // Weak edge, will be checked by hysteresis
            } else {
                edges[r][c]=0;
            }
        }
    }

    // recursive hysteresis
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (edges[r][c]==255) {
                applyHysteresis(edges, visited, r, c, lowThreshold, highThreshold);
            }
        }
    }
    
    // remove weak clusters
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (edges[r][c]!=255 && edges[r][c]!=0) {
                edges[r][c]=0;
            }
        }
    }

    // write
    ofstream outfile;
    outfile.open(outfilename);
    outfile << "P3 " << window_column << " " << window_row << " 1\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfile << edges[r][c] << " " << edges[r][c] << " " << edges[r][c] << " ";
        }
        outfile << "\n";
    }
    outfile.close();
    return 0;
}

vector<vector<int>> sobelGradientX(const vector<vector<int>>& grayscale) {
    int rows=grayscale.size();
    int cols=grayscale[0].size();
    vector<vector<int>> gradientX(rows, vector<int>(cols, 0));

    for (int y=1; y<rows-1; y++) {
        for (int x=1; x<cols-1; x++) {
            int gx=-grayscale[y-1][x-1]-2*grayscale[y][x-1]-grayscale[y+1][x-1]+grayscale[y-1][x+1]+2*grayscale[y][x+1]+grayscale[y+1][x+1];
            gradientX[y][x]=gx;
        }
    }
    return gradientX;
}

vector<vector<int>> sobelGradientY(const vector<vector<int>>& grayscale) {
    int rows=grayscale.size();
    int cols=grayscale[0].size();
    vector<vector<int>> gradientY(rows, vector<int>(cols, 0));

    for (int y=1; y<rows-1; y++) {
        for (int x=1; x<cols-1; x++) {
            int gy=-grayscale[y-1][x-1]-2*grayscale[y-1][x]-grayscale[y-1][x+1]+grayscale[y+1][x-1]+2*grayscale[y+1][x]+grayscale[y+1][x+1];
            gradientY[y][x] = gy;
        }
    }
    return gradientY;
}

vector<vector<double>> sobelDirection(vector<vector<int>>& grayscale) {
    vector<vector<int>> sobelX=sobelGradientX(grayscale);
    vector<vector<int>> sobelY=sobelGradientY(grayscale);
    vector<vector<double>> sobelDirections(window_row, vector<double>(window_column,0.0));
    
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            sobelDirections[r][c]=atan2(sobelY[r][c], sobelX[r][c]);
        }
    }
    
    return sobelDirections;
}

int part3(int argc, char *argv[]) {
    string p3or6;
    int rgb;
    int lowThreshold=180, highThreshold=190; // values 180/190 optimized for adt.jpg
    string infilename=("image.ppm"), outfilename=("image1.ppm"), outfilename2=("image2.ppm"), outfilename3=("imagef.ppm");
    
    if (argc>1) {
        for (int i=0; i+2<argc; i+=2) {
            string curr=argv[i+1];
            if (curr=="-f") {infilename=argv[i+2];}
            if (curr=="-lt") {
                string lowThreshStr=argv[i+2];
                lowThreshold=stoi(lowThreshStr);
            }
            if (curr=="-lt") {
                string highThreshStr=argv[i+2];
                highThreshold=stoi(highThreshStr);
            }
            if (curr=="-of") {outfilename=argv[i+2];}
            if (curr=="-f2") {outfilename2=argv[i+2];}
            if (curr=="-ff") {outfilename3=argv[i+2];}
        }
    }
    
    ifstream imagefile(infilename);
    imagefile >> p3or6 >> window_column >> window_row >> rgb;
    
    // GRAYSCALE
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
    
    // DOUBLE HYSTERISIS METHODS ONLY
    // sobel
    vector<vector<int>> sobel(window_row, vector<int>(window_column));
    for (int r=1; r<window_row-1; r++) {
        for (int c=1; c<window_column-1; c++) {
            int s_grad=round(sqrt(pow(grayscale[r-1][c+1]-grayscale[r+1][c-1] + 2*grayscale[r][c+1]-2*grayscale[r][c-1] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)+
                                  pow(grayscale[r+1][c-1]-grayscale[r-1][c+1] + 2*grayscale[r+1][c]-2*grayscale[r-1][c] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)));
            sobel[r][c] = s_grad;
        }
    }

    // double thresh
    vector<vector<int>> visited(window_row, vector<int>(window_column, 0));
    vector<vector<int>> edges(window_row, vector<int>(window_column, 0));
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (sobel[r][c]>highThreshold) {
                edges[r][c]=255; // Strong edge
            } else if (sobel[r][c]>lowThreshold) {
                edges[r][c]=sobel[r][c]; // Weak edge, will be checked by hysteresis
            } else {
                edges[r][c]=0;
            }
        }
    }

    // recursive hysteresis
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (edges[r][c]==255) {
                applyHysteresis(edges, visited, r, c, lowThreshold, highThreshold);
            }
        }
    }
    
    // remove weak clusters
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            if (edges[r][c]!=255 && edges[r][c]!=0) {
                edges[r][c]=0;
            }
        }
    }

    // write double thresh only
    ofstream outfile1;
    outfile1.open(outfilename);
    outfile1 << "P3 " << window_column << " " << window_row << " 1\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfile1 << edges[r][c] << " " << edges[r][c] << " " << edges[r][c] << " ";
        }
        outfile1 << "\n";
    }
    outfile1.close();
    
    // NON MAXIMUM SUPPRESSION METHODS ONLY
//     for (int r=1; r<window_row-1; r++) {
//         for (int c=1; c<window_column-1; c++) {
//             int s_grad=round(sqrt(pow(grayscale[r-1][c+1]-grayscale[r+1][c-1] + 2*grayscale[r][c+1]-2*grayscale[r][c-1] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)+
//                                   pow(grayscale[r+1][c-1]-grayscale[r-1][c+1] + 2*grayscale[r+1][c]-2*grayscale[r-1][c] + grayscale[r+1][c+1]-grayscale[r-1][c-1], 2)));
//             sobel[r][c] = s_grad;
//         }
//     }
    vector<vector<double>> direction=sobelDirection(grayscale);
    vector<vector<int>> nms(window_row, vector<int>(window_column, 0)); 
    for (int r=1; r<window_row-1; r++) {
        for (int c=1; c<window_column-1; c++) {
            double angle=direction[r][c];
            int grad=sobel[r][c];
            int grad1=0, grad2=0;

            if ((angle>=0 && angle<M_PI/8) || (angle>=15*M_PI/8 && angle<=2*M_PI) || (angle>=7*M_PI/8 && angle<9*M_PI/8)) {
                grad1=sobel[r][c+1];
                grad2=sobel[r][c-1];
            } else if ((angle>=M_PI/8 && angle<3*M_PI/8) || (angle>=9*M_PI/8 && angle<11*M_PI/8)) {
                grad1=sobel[r-1][c+1];
                grad2=sobel[r+1][c-1];
            } else if ((angle>=3*M_PI/8 && angle<5*M_PI/8) || (angle>=11*M_PI/8 && angle<13*M_PI/8)) {
                grad1=sobel[r-1][c];
                grad2=sobel[r+1][c];
            } else if ((angle>=5*M_PI/8 && angle<7*M_PI/8) || (angle>=13*M_PI/8 && angle<15*M_PI/8)) {
                grad1=sobel[r-1][c-1];
                grad2=sobel[r+1][c+1];
            }

            // Check if the current pixel is a local maximum
            if (grad>=grad1 && grad>=grad2) {
                nms[r][c]=255;
            } else {
                nms[r][c]=0;
            }
        }
    }
    
    // write nms alg only
    ofstream outfile2;
    outfile2.open(outfilename2);
    outfile2 << "P3 " << window_column << " " << window_row << " 1\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfile2 << nms[r][c] << " " << nms[r][c] << " " << nms[r][c] << " ";
        }
        outfile2 << "\n";
    }
    outfile2.close();
    
    // DOUBLE HYSTERISIS AND NMS
    vector<vector<int>> bothofthem(window_row, vector<int>(window_column, 0)); 
    for (int r=1; r<window_row-1; r++) {
        for (int c=1; c<window_column-1; c++) {
            if ((nms[r][c]==edges[r][c]) && (nms[r][c]==255)) { bothofthem[r][c]=255; }
        }
    }
    
    // write with both
    ofstream outfilef;
    outfilef.open(outfilename3);
    outfilef << "P3 " << window_column << " " << window_row << " 1\n";
    for (int r=0; r<window_row; r++) {
        for (int c=0; c<window_column; c++) {
            outfilef << bothofthem[r][c] << " " << bothofthem[r][c] << " " << bothofthem[r][c] << " ";
        }
        outfilef << "\n";
    }
    outfilef.close();
    return 0;
}

int main(int argc, char *argv[]) {
    part3(argc, argv);
}