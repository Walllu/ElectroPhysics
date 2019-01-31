#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class Prometheus {

private:
  float xmin, xmax, ymin, ymax; //Grid specifications
  unsigned N; //Resolution constant
  vector<float> potential;
  vector<float> charge;
  vector<bool> fixedpotential; //Bool type to check if potential fixed Y/N... skip on basis of condition
public:
  Prometheus(xmin, xmax, ymin, ymax, N);
  float xmin() const;
  float xmax() const;
  float ymin() const;
  float ymax() const;
  unsigned u() const;
  float h() const; //Stepsize
  
  
  
  


};



int main()
{
    cout << "Hello world!" << endl;
    int N=301; // make it odd
    double a = 5; // small radius
    double b = 10; // large radius
    double V = 100; // potential at outer cylinder

    double h = 2*b/(double)(N-1); // grid spacing
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
    }

    return 0;
}
