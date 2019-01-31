#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int dimension = 41; // x,y are going to go from [0, dimension]

struct stencil2D {
    int x_dim;
    int y_dim;
    // positve_s encodes the boundary condition locations and values (non-zero)
    double positive_s[41][41];
    // negative_s encodes the (zero) boundary condition locations (using 1's and 0's, former standing for non-zero, latter standing for the zero locations)
    double negative_s[41][41];
};

void Heat2D_next(double (*A)[41][41], double (*B)[41][41], stencil2D *stenc, bool A1st, double alpha, double dx, double dt){
    double pos[][41] = (*stenc).positive_s;
    double neg[][41] = (*stenc).negative_s;
    int j = 0;
    int k = 0;
    if(A1st){
        // A is the next, B is the previous
        for(j=0;j<dimension;j++){
            for(k=0;k<dimension;k++){
                // check stencils first
                if(pos[j][k]!=0){
                    (*A)[j][k]=pos[j][k];
                } else if(neg[j][k]==0){
                    (*A)[j][k]=pos[j][k];
                } else {
                // FTCS on A using B as previous moment
                u_j1_k_n = (*B)[j+1][k]; // this is the "j+1"th  point of the previous function
                u_j_k1_n = (*B)[j][k+1];
                u_j_1_k_n = (*B)[j-1][k];
                u_j_k_1_n = (*B)[j][k-1];
                u_j_k_n = (*B)[j][k]; // previous time step, same location
                // now for the next time step
                u_j_k_n1 = u_j_k_n + alpha*(dt/(dx*dx))*(u_j1_k_n + u_j_k1_n + u_j_1_k_n + u_j_k_1_n - 4*u_j_k_n);
                (*A)[j][k] = u_j_k_n1;
                }
            }
        }
    } else {
        // B is the next, A is the previous
        for{int j=0;j<dimension;j++}{
            for(int k=0;k<dimension;k++){
                // check stencils first
                if(pos[j][k]!=0){
                    (*B)[j][k]=pos[j][k];
                } else if(neg[j][k]==0){
                    (*B)[j][k]=pos[j][k];
                }
                // FTCS on B using A as previous moment
                u_j1_k_n = (*A)[j+1][k]; // this is the "j+1"th  point of the previous function
                u_j_k1_n = (*A)[j][k+1];
                u_j_1_k_n = (*A)[j-1][k];
                u_j_k_1_n = (*A)[j][k-1];
                u_j_k_n = (*A)[j][k]; // previous time step, same location
                // now for the next time step
                u_j_k_n1 = u_j_k_n + alpha*(dt/(dx*dx))*(u_j1_k_n + u_j_k1_n + u_j_1_k_n + u_j_k_1_n - 4*u_j_k_n);
                (*B)[j][k] = u_j_k_n1;
            }
        }
    }
}

void stencil_maker(stencil2D *stenc, int x_dim, int y_dim){
    // This function takes the multi-dimensional array and creates a stencil out of it
    double pos_arr[41][41];
    double neg_arr[41][41];
    // For now let's just make the hot walls boundary condition
    double temp = 20;
    int j = 0;
    int k = 0;
    for(j=0;j<dimension;j++){
        for(k=0;k<dimension;k++){
            // boundaries
            if(j==0 || k==0 || j==dimension-1 || k==dimension-1){
                pos_arr[j][k]=temp;
                neg_arr[j][k]=1;
            } else {
                pos_arr[j][k]=0;
                neg_arr[j][k]=0;
            }
        }
    }
    // now set the parameters in the stencil
    (*stenc).x_dim = x_dim;
    (*stenc).y_dim = y_dim;
    (*stenc).positive_s = pos_arr;
    (*stenc).negative_s = neg_arr;
}

int main(){
    // Program for calculating the FTCS Heat Diffusion from initial conditions
    double dx = 0.2;
    double dt = 0.01;
    double tmax = 10.0;
    double alpha = 1.0;
    // create stencil for simulation
    stencil2D stencil;
    stencil_maker(&stencil, dimension, dimension); // populate it with real date (hot walls, for now)
    // Create multidimensional arrays to work with
    double A[dimension][dimension];
    double B[dimension][dimension];
    bool A1st = true; // we'll switch A and B out as the "previous moment in time" and just overwrite the other
    // now let's loop through time
    double time;
    for(time=0.0;time<tmax;t+=dt){
        cout << time << endl;
        if(time==0.0){
            A = stencil;
            A1st = false; // setting it up so next time, we take A to be the previous moment in the Heat2D function
        }
        else{
            // compute the next moment with FTCS
            Heat2D_next( &A, &B, &stencil, A1st, alpha, dx, dt);
            A1st = !A1st;
        }
    }
    cout << "done!" << endl;
    if(A1st){
        cout << "B is the final result" << endl;
        cout << B[5][5] << endl;
    } else {
        cout << "A is the final result" << endl;
        cout << A[5][5] << endl;
    }
    return 0;
}