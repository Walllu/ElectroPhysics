#include <iostream>
#include <vector>
#include <cmath>
#include <string>
using namespace std;

// A simple testing of the class syntax in C++ - Walter
class Grid{
    private:
        int x_dim, y_dim; // dimensions of the grid
        //float N; // resolution
        vector<vector <float>> grid;
        float dx, dy;
    public:
        Grid() : x_dim(11), y_dim(11), dx(0.2), dy(0.2) {
            cout << "hellloooo?" << endl;
            vector<float> temp(x_dim);
            int i;
            vector< vector<float> > temp2;
            for(i=0;i<x_dim;i++){
                temp2.push_back(temp);
            }
            grid = temp2;
            //vector< vector<float> > temp3(x_dim, 0.0);
            //grid = temp3;
        }
        Grid(int num): x_dim(41), y_dim(41), dx(0.2), dy(0.2){
            Grid grid; // make default grid
            // now we want to create problem 1
            cout << "hi hi" << endl;
            if(num==1){
                cout << "this is grid 1" << endl;
            }
        }
        void get(){
            cout << x_dim << endl;
            for(int i=0; i<x_dim; i++){
                for(int j=0; j<y_dim; j++){
                  //  cout << grid[i][j];
                }
                //cout << endl;
            }
            cout << grid.size() << endl;
        }
};




int main(){
    // testing the Grid class
    cout << "hey?" << endl;
    Grid grid;
    grid.get();
    Grid grid2(1);    
    return 0;
}