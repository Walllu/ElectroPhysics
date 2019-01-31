#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// A simple testing of the class syntax in C++ - Walter
class Grid{
    private:
        int x_dim, y_dim; // dimensions of the grid
        float N; // resolution
        vector<vector <float>> grid;
        float dx, dy;
    public:
        Grid() {
            x_dim = 40;
            y_dim = 40;
            N = 10; // represents the physical dimensions of the grid
            dx = N/x_dim;
            dy = N/y_dim;
            vector<float> temp(x_dim);
            int i;
            vector< vector<float> > temp2;
            for(i=0;i<x_dim;i++){
                temp2.push_back(temp);
            }
            grid = temp2;
        }
        void get(){
            cout << x_dim << endl;
            //cout << grid << endl;
        }
};



int main(){
    // testing the Grid class
    Grid grid;
    grid.get();
    return 0;
}