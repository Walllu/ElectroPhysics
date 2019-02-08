#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// a helper class that defines the grid and grid spacing
class Grid{
    private:
        int x_dim, y_dim; //  number of nodes in x, y direction, (x, y index goes from 0 to x_dim - 1, y_dim - 1)
        float dx, dy; // grid spacing in each x,y direction
        float x_dist, y_dist; // physical x and y distance
        vector<vector <double> > potential; // the actual grid of some dimension
        vector<vector <bool> > constants;
    public:
        /*Grid() : x_dim(11), y_dim(11), dx(0.2), dy(0.2) {
            // default configuration
            vector<double> y_axis(y_dim); // construct the grid of zeros
            vector< vector<double> > x_axis;
            for(int i=0;i<x_dim;i++){
                x_axis.push_back(y_axis);
            }
            grid = x_axis; // default 0's matrix
        }*/
        Grid(int x_dim_arg, int y_dim_arg, float x_dist_arg, float y_dist_arg): x_dim(x_dim_arg), y_dim(y_dim_arg), x_dist(x_dist_arg), y_dist(y_dist_arg){
            // this is a more flexible default starter grid of zeros
            // just define the x,y array size and physical depth of each
            // and we set dx, dy here based on physical_length/grid_spaces
            dx = x_dist/(x_dim - 1);
            dy = y_dist/(y_dim - 1);

            vector <double> y_axis_1(y_dim); // construct the grid of zeros
            vector < vector<double> > x_axis_1;
            for(int i=0; i<x_dim; i++){
                x_axis_1.push_back(y_axis_1);
            }
            potential = x_axis_1; // default 0's matrix

            vector <bool> y_axis_2(y_dim); // construct the grid of zeros
            vector < vector<bool> > x_axis_2;
            for(int i=0; i<x_dim; i++){
                x_axis_2.push_back(y_axis_2);
            }
            constants = x_axis_2;
        }
        int get_x_dim(){
            return x_dim; // getter for x dimension
        }
        int get_y_dim(){
            return y_dim; // getter for y dimension
        }
        float get_dx(){
            return dx; // getter for dx space-step
        }
        float get_dy(){
            return dy; // getter for dy space-step
        }
        float get_x_dist(){
            return x_dist; // getter for dy space-step
        }
        float get_y_dist(){
            return y_dist; // getter for dy space-step
        }
        float get_x_position(int i){
            return i*dx;
        }
        float get_y_position(int i){
            return i*dy;
        }
        void set_point(int i,int j, double val){
            potential[i][j] = val; // sets a point in the grid to "val"
            // i refers to the y position, j to the x position
        }
        void set_constant(int i,int j){
            constants[i][j] = true; // sets a point in the grid to "val"
            // i refers to the y position, j to the x position
        }
        double get_point(int i,int j){
            return potential[i][j]; // getter for value of point in grid
            // i refers to the y position, j to the x position
        }
        bool is_it_constant(int i, int j){
            return constants[i][j];
        }
        void grid_print(){  // method for printing the whole grid
            for(int i=0; i<x_dim; i++){
                for(int j=0; j<y_dim; j++){
                    cout << "(" << potential[i][j] << ")";
                }
                cout << endl;
            }
        }
};
// function for setting initial conditions
void Set_initial (int num, Grid& grid){
        if(num==1){
                //create problem 1, I'll just do hot walls)
                for(int i=0; i<grid.get_x_dim(); i++){
                    for(int j=0; j<grid.get_y_dim(); j++){
                        if(i==0||j==0||i==grid.get_y_dim()-1||j==grid.get_x_dim()-1){
                            //A[i][j]=20; // location and value of point
                            //B[i][j]=1; // mark as constant
                            grid.set_point(i,j,20);
                            grid.set_constant(i,j);
                        } else{
                            //A[i][j]=0;  //location and value of point
                            //B[i][j]=0; // mark as variable
                            grid.set_constant(i,j);
                        }
                    }
                }
            }
            else if(num==2){
                for(int i=0; i<grid.get_x_dim(); i++){
                    for(int j=0; j<grid.get_y_dim(); j++){
                        if(j==0){
                            //A[i][j]=20; // location and value of point
                            //B[i][j]=1; // mark as constant
                            grid.set_point(i,j,20);
                            grid.set_constant(i,j);
                        } else if(j==grid.get_y_dim()-1){
                            grid.set_point(i,j,-20);
                            grid.set_constant(i,j);
                        }
                        else{
                            //A[i][j]=0;  //location and value of point
                            //B[i][j]=0; // mark as variable
                            grid.set_point(i,j,0);
                            //constant_locations.set_point(i,j,0);
                        }
                    }
                }
            }
        }

// function for evaluating the next time-step via FTCS
void Heat2D_next_u(Grid& next_u, Grid& prev_u, const double alpha, float dt){
    // Now we want to loop through all the spatial coordinates of the next_u, contructing it
    int x_ext = next_u.get_x_dim();
    int y_ext = next_u.get_y_dim();
    float dx = next_u.get_dx();
    float dy = next_u.get_dy(); // setting up all the values we need
    double courant_factor1 = 2.0*alpha*(dt/pow(dx,2));
    double courant_factor2 = 2.0*alpha*(dt/pow(dy,2));
    for(int j=0;j<x_ext;j++){
        for(int k=0;k<y_ext;k++){ // gonna switch over to "k" for my convenience - it's just an index
            // check if the point is constant
            if(prev_u.is_it_constant(j,k)){
                // if it is constant, we want to ignore it.
                continue; // Move onto next iteration
                //next_u.set_point(j,k,stencil.get_value_at_point(j,k)); // sets the value as constant value from stencil
            } else if(j==0||k==0||j==x_ext-1||k==y_ext-1){
                // if we have an edge case that is not constant
                double u_j1_k_n,  u_j_k1_n, u_j_1_k_n, u_j_k_1_n, u_j_k_n, u_j_k_n1;
                int k_1, j_1, k1, j1;
                k_1 = k - 1;
                j_1 = j - 1;
                k1 = k + 1;
                j1 = j + 1;
                // if j=0, then a decrease in k is bad
                if(j==0){
                    j_1 = x_ext - 1;
                }
                // if k=0, then a decrease in j is bad
                if(k==0){
                    k_1 = y_ext - 1;
                }
                // if j==x_ext-1 then an increase in k is bad
                if(j==x_ext-1){
                    j1 = 0;
                }
                // if k--y_ext-1 then an increase in j is bad
                if(k==y_ext-1){
                    k1 = 0;
                }
                u_j1_k_n = prev_u.get_point(j1, k);
                u_j_k1_n = prev_u.get_point(j, k1);
                u_j_1_k_n = prev_u.get_point(j_1,k);
                u_j_k_1_n = prev_u.get_point(j,k_1);
                u_j_k_n = prev_u.get_point(j,k); // now we have all the points from prev_u
                u_j_k_n1 = u_j_k_n + courant_factor1*(u_j1_k_n -2*u_j_k_n + u_j_1_k_n) + courant_factor2*(u_j_k_1_n + u_j_k1_n - 2*u_j_k_n);
                next_u.set_point(j,k, u_j_k_n1);
            } else { // if not constant, then we run FTCS on that point
                double u_j1_k_n,  u_j_k1_n, u_j_1_k_n, u_j_k_1_n, u_j_k_n, u_j_k_n1;
                // as variable names can't contain '+', '-', or '[]' characters
                // I've made the following notation:
                    // if u_{j,k}^n denotes the (j,k) point at time-step 'n'
                    // --> u_{j+1,k}^n  is  u_j1_k_n
                    // --> u_{j,k-1}^n  is  u_j_k_1_n
                u_j1_k_n = prev_u.get_point(j+1, k);
                u_j_k1_n = prev_u.get_point(j, k+1);
                u_j_1_k_n = prev_u.get_point(j-1,k);
                u_j_k_1_n = prev_u.get_point(j,k-1);
                u_j_k_n = prev_u.get_point(j,k); // now we have all the points from prev_u
                u_j_k_n1 = u_j_k_n + courant_factor1*(u_j1_k_n -2*u_j_k_n + u_j_1_k_n) + courant_factor2*(u_j_k1_n + u_j_k_1_n - 2*u_j_k_n);
                next_u.set_point(j,k, u_j_k_n1);
            }
        }
    }
}
void timeloop(){

    Grid A(9, 11, 1, 1);
    // these two are going to be switiching back and forth
    Set_initial(2, A);
    Grid B = A;

    bool A_next = true; // this boolean is going to determine whether A or B is going to store the next moment in time
    // working with default initial conditions dx, dy, dt, tmax
    int maximum_time_in_seconds, transform;
    float dt;
    // time is an integer, time step-size (dt) is a float, and transform is an integer used for debugging purposes - enter 0 for transform if you don'e give a shit about debugging
    cout << "Maximum time in seconds, and time step-size (dt) : transform too" << endl;
    cin >> maximum_time_in_seconds >> dt >> transform;  // "transform" is the iteration number that you want to view - transform = 1 is the first time step, 2 is 2nd, etc.
    cout << "Input your value for alpha" << endl;
    float alpha;
    cin >> alpha;
    //cout << " Courant Condition?" << alpha*(dt/pow(A.get_dx(), 2)) << endl;
    cout << "Starting FTCS with tmax: " << maximum_time_in_seconds << ", dt: " << dt << ", and alpha: " << alpha << endl;
    // Now let's begin looping through time
    for(double time = 0.0; time<maximum_time_in_seconds; time+=dt){
        if(time==0.0){
            A_next = !A_next; // flip the boolean
        } else {
            if(A_next){
                // A is next up --> next_u = A; prev_u = B;
                Heat2D_next_u(A, B, alpha, dt);
            } else{
                // B is next up --> next_u = B; prev_u = A;
                Heat2D_next_u(B, A, alpha, dt);
            }
            if(time == transform*dt){ // debugging purposes
                cout << time << endl;
                A.grid_print();
                cout <<"A^^"<< endl;
                B.grid_print();
                cout <<"B^^"<< endl;
            }
            A_next = !A_next;
        }
    }

    cout << "Done! : " << A_next << endl;
    if(A_next){
        B.grid_print();
    } else {
         A.grid_print();
    }
}
int main(){
    // Program for calculating the FTCS Heat Diffusion from initial conditions - modifying the initial conditions and calling the correct Stencil constructor
    // one can recreate an electric potential situation as well.
    // The program is structured so that we can run a timeloop on a specfic situation by creating the Stencil, and calling the "timeloop" function on it
    // for example:

    // create stencil for the simulation
     //stenc.stencil_print();
    // call the timeloop
    //stenc.stencil_print();
    timeloop();

    // repeat the above for other initial conditions.
    return 0;
}
