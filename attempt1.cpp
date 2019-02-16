#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class Prometheus {

private:
  double xmin, xmax, ymin, ymax; //Grid specifications
  // as for now, keep the x and y distances the same
  unsigned N; // number of nodes in each direction, ie index goes from 0 to N-1
  vector<float> potential;
  vector<float> charge;
  vector<bool> fixedpotential; //Bool type to check if potential fixed Y/N... skip on basis of condition
public:
  Prometheus(double x1, double x2, double y1, double y2, unsigned index_scope);
  void set_initial();
  //double xmin() const;
  //double xmax() const;
  //double ymin() const;
  //double ymax() const;
  //unsigned u() const;
  double h() const; //Stepsize
  double x(unsigned x_index);
  double y(unsigned y_index);

};

// Constructor to read grid specifications
Prometheus::Prometheus(double x1, double x2, double y1, double y2, unsigned index_scope) {
    xmin = x1;
    xmax = x2;
    ymin = y1;
    ymax = y2;
    N = index_scope;
};

// h returns grid spacing
double Prometheus::h() const {
    return (xmax-xmin)/(N-1);
}

// x returns horizontal position as function of horizontal index
double Prometheus::x(unsigned x_index){
    return xmin + h()*x_index;
};

// y returns vertical position as a function of vertical index
double Prometheus::y(unsigned y_index){
    return ymin + h()*y_index;
};

// This function is supposed to set a initial grid
void Prometheus::set_initial(){
    for (unsigned i = 0; i<N; i++){
            for (unsigned i = 0; i<N; i++){
            }
    }
};



int main()
{
    cout << "Hello world!" << endl;

    Prometheus simulation(0, 10, 0, 10, 101);
    simulation.set_initial();




    /*double h = 2*b/(double)(N-1); // grid spacing
    // i - vertical index, j - horizontal index, 0<=i,j<N
    int p; // p = i + jN


    vector<double> grid(N*N, 0);
    vector< vector<double> > iterations;

    int di=0;
    int dj=0;

    vector<int> out;
    vector<double> next;

    //for (p = 0; p < N*N; p++) {
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){

            di = i - 10;
            dj = j - 10;
            if (abs(a - h*sqrt(di*di + dj*dj)) <= h/2.0) {
                grid[i+j*N] = 0;
            }
            if (abs(b - h*sqrt(di*di + dj*dj)) <= h/2.0) {
                grid[i+j*N] = V;
            }
            if (h*sqrt(di*di + dj*dj) >= b - h/2.0){
                out.push_back(i+j*N);
            }
            if (h*sqrt(di*di + dj*dj) <= a + h/2.0){
                out.push_back(i+j*N);
            }
        }
    }


    iterations.push_back(grid);

    double rel = 1;
    if (rel > 0.01){
        for (p=0; p<N*N; p++){

            if(find(out.begin(), out.end(), p) != out.end()) {
                next.push_back(iterations.back()[p]);
            }
            else{
                next.push_back(iterations.back()[p] + (1/2)*(iterations.back()[p+1] + iterations.back()[p-1] + iterations.back()[p+N]));
                if (abs(next.back()-iterations.back()[p])/iterations.back()[p] > rel){
                    rel = abs(next.back()-iterations.back()[p])/iterations.back()[p];
                }
            }

        }
        iterations.push_back(next);
        next.clear();
    }*/

    return 0;
}
