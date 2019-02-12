#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// a helper class that defines the grid and grid spacing
class Grid{
    private:
        int x_dim, y_dim; //  number of nodes in x, y direction, (x, y index goes from 0 to x_dim - 1)
        float x_phys, y_phys; //
        float dx, dy;   // grid spacing in each x,y direction
        vector<vector <double> > grid; // the actual grid of some dimension
    public:
        Grid() : x_dim(13), y_dim(13), x_phys(1), y_phys(1) {
            // default configuration
            dx = x_phys/x_dim;
            dy = y_phys/y_dim;
            vector<double> y_axis(y_dim); // construct the grid of zeros
            vector< vector<double> > x_axis;
            for(int i=0;i<x_dim;i++){
                x_axis.push_back(y_axis);
            }
            grid = x_axis; // default 0's matrix
        }
        Grid(int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg): x_dim(x_dim_arg), y_dim(y_dim_arg), x_phys(x_phys_arg), y_phys(y_phys_arg) {
            // this is a more flexible default starter grid of zeros
            // just define the x,y array size and physical depth of each
            // and we set dx, dy here based on physical_length/grid_spaces
            dx = x_phys/(x_dim-1);
            dy = y_phys/(y_dim-1);
            vector <double> y_axis(y_dim); // construct the grid of zeros
            vector < vector<double> > x_axis;
            for(int i= 0; i<x_dim; i++){
                x_axis.push_back(y_axis);
            }
            grid = x_axis; // default 0's matrix
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
        float get_x_position(int i){
            return i*dx - x_phys/2.0;
        }
        float get_y_position(int i){
            return i*dy - y_phys/2.0;
        }
        void set_point(int i,int j, double val){
            grid[i][j] = val; // sets a point in the grid to "val"
            // i refers to the y position, j to the x position
        }
        double get_point(int i,int j){
            return grid[i][j]; // getter for value of point in grid
            // i refers to the y position, j to the x position
        }
        void grid_print(){  // method for printing the whole grid
            for(int i= 0; i<x_dim; i++){
                for(int j= 0; j<y_dim; j++){
                    cout << "(" << grid[i][j] << ")";
                }
                cout << endl;
            }
        }
};
// class for making the initial value stencils for problems
class Stencil{
    private:
        // Grids that represnet the locations and values of all the constants
        Grid problem_values, problem_constant_locations;
        // problem_values are the 'locations and values of all points in the initial moment
        // problem_constant_locations has the locations of all constant values in the grid
    public:
        Stencil(int num){ // we make a constructor that acts differently depending on which input you start with
            if(num==1){
                int p = 13;
                int r = 13;
                float m = 1.0;
                float n = 1.0;
                //create problem 1, I'll just do hot walls)
                Grid initial_values(p, r, m, n), constant_locations(p, r, m, n);
                for(int i = 0; i < initial_values.get_x_dim(); i++){
                    for(int j= 0 ; j < initial_values.get_y_dim(); j++){
                        if(i==initial_values.get_x_dim()-1||j==initial_values.get_y_dim()-1||i== 0||j== 0){

                            initial_values.set_point(i,j,20);
                            constant_locations.set_point(i,j,1);
                        } else{

                            initial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,0);
                        }
                    }
                }
                problem_values = initial_values;
                problem_constant_locations = constant_locations; // 1 represents constant, 0 represents variable
                // end of problem 1
            }
            else if(num==2){
                int p = 13;
                int r = 13;
                float m = 1.0;
                float n = 1.0;
                //create problem 1, I'll just do hot walls)
                Grid intitial_values(p, r, m, n), constant_locations(p, r, m, n);
                for(int i = 0; i < intitial_values.get_x_dim(); i++){
                    for(int j=0; j < intitial_values.get_y_dim(); j++){
                        if(j==intitial_values.get_y_dim()-1||j== 0){

                            intitial_values.set_point(i,j,20);
                            constant_locations.set_point(i,j,1);
                        } else{

                            intitial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,0);
                        }
                    }
                }
                problem_values = intitial_values;
                problem_constant_locations = constant_locations; // 1 represents constant, 0 represents variable
                // end of problem 1
            }
            else if(num==3){
                int p = 13;
                int r = 13;
                float m = 1.0;
                float n = 1.0;
                //create problem 1, I'll just do hot walls)
                Grid initial_values(p, r, m, n), constant_locations(p, r, m, n);
                float x, y;
                float a, b;
                a = 0.5;
                b = 0.2;
                for(int i = 0; i < initial_values.get_x_dim(); i++){
                    x = initial_values.get_x_position(i);
                    for(int j= 0; j < initial_values.get_y_dim(); j++){
                        y = initial_values.get_y_position(j);
                        if(x*x + y*y >= a*a){

                            initial_values.set_point(i,j,20);
                            constant_locations.set_point(i,j,1);
                        }
                        else if(x*x + y*y <= b*b){

                            initial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,1);
                        }
                        else{

                            initial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,0);
                        }
                    }
                }
                problem_values = initial_values;
                problem_constant_locations = constant_locations; // 1 represents constant, 0 represents variable
                // end of problem 1
            }
            else if(num==4){
                int p = 13;
                int r = 13;
                float m = 1.0;
                float n = 1.0;
                //create problem 1, I'll just do hot walls)
                Grid initial_values(p, r, m, n), constant_locations(p, r, m, n);
                float x, y;
                float a, b;
                a = 1;
                b = 0.2;
                for(int i = 0; i < initial_values.get_x_dim()-1; i++){
                    x = initial_values.get_x_position(i);
                    for(int j= 0; j < initial_values.get_y_dim(); j++){
                        y = initial_values.get_y_position(j);
                        if(j = 0){

                            initial_values.set_point(i,j,20);
                            constant_locations.set_point(i,j,1);
                        }
                        if(j = initial_values.get_y_dim()-1){

                            initial_values.set_point(i,j,-20);
                            constant_locations.set_point(i,j,1);
                        }
                        else if(x*x + y*y <= b*b){

                            initial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,1);
                        }
                        else{

                            initial_values.set_point(i,j,0);
                            constant_locations.set_point(i,j,0);
                        }
                    }
                }
                problem_values = initial_values;
                problem_constant_locations = constant_locations; // 1 represents constant, 0 represents variable
                // end of problem 1
            }




        }
        void stencil_print(){
            for(int i=0;i<problem_values.get_x_dim(); i++){
                for(int j=0;j<problem_values.get_y_dim(); j++){
                    cout << "(" << problem_values.get_point(i,j) << "," << problem_constant_locations.get_point(i,j) << ")";
                }
                cout << endl;
            }
        }
        Grid get_values(){
            return problem_values;
        }
        Grid get_constants(){
            return problem_constant_locations;
        }
        double get_value_at_point(int i, int j){
            return problem_values.get_point(i,j);
        }
        bool is_it_constant(int i, int j){
            if(problem_constant_locations.get_point(i,j)==1){
                return true;
            } else {
                return false;
            }
        }
};
// function for evaluating the next time-step via FTCS
void Heat2D_next_u(Grid& next_u, Grid& prev_u, Stencil& stencil, const double alpha, float dt){
    // Now we want to loop through all the spatial coordinates of the next_u, contructing it
    int x_ext = next_u.get_x_dim();
    int y_ext = next_u.get_y_dim();
    float dx = next_u.get_dx();
    float dy = next_u.get_dy(); // setting up all the values we need
    double courant_factor1 = 2.0*alpha*(dt/pow(dx,2));
    double courant_factor2 = 2.0*alpha*(dt/pow(dy,2));
    for(int j=0;j<x_ext;j++){
        for(int k= 0;k<y_ext;k++){ // gonna switch over to "k" for my conventience - it's just an index
            // check if the point is constant
            if(stencil.is_it_constant(j,k)){
                // if it is constant, we want to ignore it.
                continue; // Move onto next iteration
                //next_u.set_point(j,k,stencil.get_value_at_point(j,k)); // sets the value as constant value from stencil
            } else if(j==0||k== 0||j==x_ext-1||k==y_ext-1){
                // if we have an edge case that is not constant
                double u_j1_k_n,  u_j_k1_n, u_j_1_k_n, u_j_k_1_n, u_j_k_n, u_j_k_n1;
                int k_1, j_1, k1, j1;
                k_1 = k - 1;
                j_1 = j - 1;
                k1 = k + 1;
                j1 = j + 1;
                // if j=0, then a decrease in k is bad
                if(j==0){
                    j_1 = x_ext-1;
                }
                // if k=0, then a decrease in j is bad
                if(k == 0){
                    k_1 = y_ext-1;
                }
                // if j==x_ext-1 then an increase in k is bad
                if(j==x_ext - 1){
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
                u_j_k_n1 = u_j_k_n + courant_factor1*(u_j1_k_n -2*u_j_k_n + u_j_k1_n) + courant_factor2*(u_j_1_k_n + u_j_k_1_n - 2*u_j_k_n);
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
                u_j_k_n1 = u_j_k_n + courant_factor1*(u_j1_k_n -2*u_j_k_n + u_j_k1_n) + courant_factor2*(u_j_1_k_n + u_j_k_1_n - 2*u_j_k_n);
                next_u.set_point(j,k, u_j_k_n1);
            }
        }
    }
}
void timeloop(Stencil& stencil){

    // implements a timeloop for each initial condition stencil
    /*int p, r;
    p = 11;
    r = 11;
    float m, n;
    m = 1;
    n = 1;*/
    Grid A, B;
    // these two are going to be switiching back and forth

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
            // set the initial conditions up
            A = stencil.get_values(); // one time allocations of memeory
            B = A;
            A_next = !A_next; // flip the boolean
        } else {
            if(A_next){
                // A is next up --> next_u = A; prev_u = B;
                Heat2D_next_u(A, B, stencil, alpha, dt);
            } else{
                // B is next up --> next_u = B; prev_u = A;
                Heat2D_next_u(B, A, stencil, alpha, dt);
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
    Stencil stencil(3);
    stencil.stencil_print();
    timeloop(stencil);

    // repeat the above for other initial conditions.
    return 0;
}
