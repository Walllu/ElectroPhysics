#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// a helper class that defines the grid and grid spacing
class Grid {
	protected:
		int x_dim, y_dim; //  number of nodes in x, y direction, (x, y index goes from 0 to x_dim - 1)
		float x_phys, y_phys; // physical x, y distances
		float dx, dy;   // grid spacing in each x,y direction
		vector<vector <double> > grid; // the actual grid of some dimension
	public:
		// empty default constructor
		Grid(){}
		// constructor for making Grids defined by x,y dimension, and x,y physical dimension
		Grid(int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg) : x_dim(x_dim_arg), y_dim(y_dim_arg), x_phys(x_phys_arg), y_phys(y_phys_arg) {
			dx = x_phys / (x_dim - 1);
			dy = y_phys / (y_dim - 1);
			vector <double> y_axis(y_dim); // construct the grid of zeros
			vector < vector<double> > x_axis;
			for (int i = 0; i < x_dim; i++) {
				x_axis.push_back(y_axis);
			}
			grid = x_axis; // default 0's matrix
		}
		int get_x_dim() {
			return x_dim; // getter for x dimension
		}
		int get_y_dim() {
			return y_dim; // getter for y dimension
		}
		float get_dx() {
			return dx; // getter for dx space-step
		}
		float get_dy() {
			return dy; // getter for dy space-step
		}
		float get_x_phys() {
			return x_phys; // getter for x_phys physical dimenion
		}
		float get_y_phys() {
			return y_phys; // getter for y_phys physical dimenion
		}
		float get_x_position(int i) {
			return i * dx - x_phys / 2.0; // convert x index to x coordinate
		}
		float get_y_position(int i) {
			return i * dy - y_phys / 2.0; // convert y index to y coordinate
		}
		void set_point(int i, int j, double val) {
			grid[i][j] = val; // sets a point in the grid to "val"
			// i refers to the y position, j to the x position
		}
		double get_point(int i, int j) {
			return grid[i][j]; // getter for value of point in grid
			// i refers to the y position, j to the x position
		}
};

class MMatrix : public Grid {
	public:
		// This constructor works specifically to create an 'M' matrix specifically from a Grid
		// just call this constructor (inputting a Grid) and it will work
		MMatrix(Grid grid_to_construct_from) {
			x_dim = grid_to_construct_from.get_x_dim(); // set x_dim of 'M' matrix
			y_dim = grid_to_construct_from.get_y_dim(); // set y_dim of 'M' matrix - could be different in general
			int M_dim = x_dim * y_dim; // helper variable for generating grid
			cout << M_dim << endl;
			vector< vector<double> > Mcolumn; // now let's create this thing...
			// for each grid point in the original we want to make a new 'M' matrix row
			for (int current_point = 0; current_point < M_dim; current_point++) {
				vector<double> Mrow(M_dim);
				Mrow[current_point] = -4; // this will always occur
				// if point is in first row or column, or in the last row or column...
				// basically, let's take care of the edge cases if they occur at all
				if (current_point < x_dim) {
					// current point is in first row
					Mrow[current_point + x_dim] = 1; // below
					Mrow[M_dim-x_dim + current_point] = 1; // ------- pacman - first row
					if (current_point == 0) { // -------------------first row, first col
						Mrow[current_point + 1] = 1;
						Mrow[x_dim-1] = 1; //--------- -------------- pacman - first row
					}
					else if (current_point == x_dim - 1) { //--------first row, last col
						Mrow[current_point - 1] = 1;
						//Mrow[M_dim-1] = 1; // -------------------- pacman - first row
						Mrow[0] = 1; //----------------------------- pacman - first row
					}
					else {
						Mrow[current_point + 1] = 1;
						Mrow[current_point - 1] = 1;
					}
				}
				else if (current_point >= M_dim - x_dim) {
					// current point is in last row
					Mrow[current_point - x_dim] = 1; // above
					Mrow[current_point - M_dim + x_dim] = 1;//------- pacman - last row
					if (current_point == M_dim - x_dim) { //---------last row, first col
						Mrow[current_point + 1] = 1;
						Mrow[M_dim-1] = 1; //------------------------ pacman - last row
						//Mrow[0] = 1;//----------------------------- pacman - last row
					}
					else if (current_point == M_dim - 1) {//---------last row, last col
						Mrow[current_point - 1] = 1;
						Mrow[M_dim-x_dim] = 1; //------------------- pacman - last row
						//Mrow[x_dim-1] = 1;//---------------------- pacman - last row
					}
					else {
						Mrow[current_point + 1] = 1;
						Mrow[current_point - 1] = 1;
					}
				}
				else if (current_point % x_dim == 0 || current_point % x_dim == x_dim - 1) {
					Mrow[current_point + x_dim] = 1; // below
					Mrow[current_point - x_dim] = 1; // above
					if (current_point % x_dim == 0) {
						// do something if in first column in last row
						Mrow[current_point + x_dim - 1] = 1; //pacman
						Mrow[current_point + 1] = 1;
						Mrow[current_point + x_dim-1] = 1; //--------- pacman - first col
					}
					else {
						Mrow[current_point - 1] = 1;
						Mrow[current_point - x_dim+1] = 1; //---------- pacman - last col
					}

				}
				else {
					Mrow[current_point + x_dim] = 1; // below
					Mrow[current_point - x_dim] = 1; // above
					Mrow[current_point + 1] = 1; // right
					Mrow[current_point - 1] = 1; // left
				}
				// now that we've constructed the row (corresponding to one grid point!) we push to columns stack
				Mcolumn.push_back(Mrow);
			};
			// now we've completed the 'M' matrix, let the grid be this matrix
			grid = Mcolumn;
			x_dim = M_dim;
			y_dim = M_dim;
		}
};

// class for making the initial value stencils for problems
class Stencil {
	private:
		// Grids that represent the locations and values of all the constants
		Grid problem_values, problem_constant_locations;
		// problem_values are the 'locations and values of all points in the initial moment
		// problem_constant_locations has the locations of all constant values in the grid
	public:
		Stencil(int num) { // we make a constructor that acts differently depending on which input you start with
			cout << "Input number of x, y dimensions and x, y distances" << endl;
			int x_dim, y_dim;
			float x_phys, y_phys;
			cin >> x_dim >> y_dim >> x_phys >> y_phys;
			Grid initial_values(x_dim, y_dim, x_phys, y_phys), constant_locations(x_dim, y_dim, x_phys, y_phys);
			if (num == 1) {
			//create problem 1, I'll just do hot walls)
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					if (i == initial_values.get_x_dim() - 1 || j == initial_values.get_y_dim() - 1 || i == 0 || j == 0) {

						initial_values.set_point(i, j, 20);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
			}
			} else if (num == 2) {
			//create problem 1, I'll just do hot walls)
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					if (j == initial_values.get_y_dim() - 1) {

						initial_values.set_point(i, j, 20);
						constant_locations.set_point(i, j, 1);
					} else if(j == 0){
						initial_values.set_point(i, j, -20);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
			}
			} else if (num == 3) {
			//create problem 3, two cylinders)
			float x, y;
			float out_r, inn_r;
			cout << "type inner and outer radius" << endl;
			cin >> inn_r >> out_r;
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				x = initial_values.get_x_position(i);
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					y = initial_values.get_y_position(j);
					if (x*x + y * y >= out_r * out_r) {

						initial_values.set_point(i, j, 20);
						constant_locations.set_point(i, j, 1);
					}
					else if (x*x + y * y <= inn_r * inn_r) {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
				}
			} else if (num == 4) {
			//create problem 4, cylinder between walls)
			float x, y;
			float out_r, inn_r;
			cout << "type inner radius" << endl;
			cin >> inn_r;
			for (int i = 0; i < initial_values.get_x_dim() - 1; i++) {
				x = initial_values.get_x_position(i);
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					y = initial_values.get_y_position(j);
					if (j = 0) {

						initial_values.set_point(i, j, 20);
						constant_locations.set_point(i, j, 1);
					}
					if (j = initial_values.get_y_dim() - 1) {

						initial_values.set_point(i, j, -20);
						constant_locations.set_point(i, j, 1);
					}
					else if (x*x + y * y <= inn_r * inn_r) {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
			}
			} else if (num == 5) {
			//create problem 5, coupled edge)
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					if (j == initial_values.get_y_dim() - 1) {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 1);
					}
					else if (j == 0) {
						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
			}
			double a, b; // x y size of the rectangles
			cout << endl;
			cout << "Type x, y size of the rectangles" << endl;
			cin >> a, b;
			double V;
			cout << "provide potential" << endl;
			cin >> V;
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
			for (int j = 0; j < initial_values.get_y_dim(); j++) {
			if (abs(initial_values.get_y_position(j)) <= b / 2.0) {
				if (abs(initial_values.get_x_position(i) - initial_values.get_x_phys() / 8.0) <= a / 2.0) {
					initial_values.set_point(i, j, V);
					constant_locations.set_point(i, j, 1);
				}
				else if (abs(initial_values.get_x_position(i) + initial_values.get_x_phys() / 8.0) <= a / 2.0) {
					initial_values.set_point(i, j, -V);
					constant_locations.set_point(i, j, 1);
				}
				else if (abs(initial_values.get_x_position(i) - 3.0*initial_values.get_x_phys() / 8.0) <= a / 2.0) {
					initial_values.set_point(i, j, 0);
					constant_locations.set_point(i, j, 1);
				}
				else if (abs(initial_values.get_x_position(i) + 3.0*initial_values.get_x_phys() / 8.0) <= a / 2.0) {
					initial_values.set_point(i, j, 0);
					constant_locations.set_point(i, j, 1);
				}
			}
			}
			}

		}
			problem_values = initial_values;
			problem_constant_locations = constant_locations; // 1 represents constant, 0 represents variable
		}
		Grid get_values() {
			return problem_values;
		}
		Grid get_constants() {
			return problem_constant_locations;
		}
		double get_value_at_point(int i, int j) {
			return problem_values.get_point(i, j);
		}
		bool is_it_constant(int i, int j) {
			if (problem_constant_locations.get_point(i, j) == 1) {
				return true;
			}
			else {
				return false;
			}
		}
};

// class for linearized vector
class Linear {
	private:
		vector<double> vec_potential; //linearized potential vector
		vector<int> vec_const; //linearized vector of constant points
		int size; // length of the arrays
		int x_dim; // x dim
		int y_dim; // y dim
	public:
		Linear(Stencil s) {
			x_dim = s.get_values().get_x_dim();
			y_dim = s.get_values().get_y_dim();
			size = x_dim * y_dim;
			vector<double> temp_vec_potential(size, 0);
			for (int i = 0; i < x_dim; i++) {
				for (int j = 0; j < x_dim; j++) {
					temp_vec_potential[i + j * x_dim] = s.get_value_at_point(i, j);
				}
			}
			vec_potential = temp_vec_potential;

			vector<int> make_vec_const(size, 0);
			for (int i = 0; i < x_dim; i++) {
				for (int j = 0; j < x_dim; j++) {
					if (s.is_it_constant(i, j)) {
						make_vec_const[i + j * x_dim] = 1;
					}
					else {
						make_vec_const[i + j * x_dim] = 0;
					}
				}
			}
			vec_const = make_vec_const;
		};
		float get_size() {
			return size;
		}
		int get_x_dim(){
			return x_dim;
		}
		int get_y_dim(){
			return y_dim;
		}
		double get_value_linear(int i) {
			return vec_potential[i];
		}
		void set_value_linear(int i, double val) {
			vec_potential[i] = val;
		}
		bool is_it_constant_linear(int i) {
			if (vec_const[i] == 1) {
				return true;
			}
			else {
				return false;
			}
		}
};

// iteration function for obtain next Linear class object from previous one
bool next_Jacobi(Linear& prev, Linear& next, MMatrix& mat, double tolerate) {
	double temp; // a temporary variable that stores the result of matrix multiplication
	int size = prev.get_size();
	int x_dim = prev.get_x_dim();
	int y_dim = prev.get_y_dim();
	double largest_change = 0; // store the largest absolute change for a grid point in iteration
	bool keep_iterating = true; // this boolean will inform the time loop to keep iterating, given tolerance
	for (int n = 0; n < size; n++) {
		// if the current point is not constant...
		if (!(prev.is_it_constant_linear(n))) {
			temp = 0; // need to reset temp value to zero
			for (int m = 0; m < size; m++) {
				// loop over matrix to do multiplication
				if (m != n) {
					temp += prev.get_value_linear(m)*mat.get_point(n, m);
				}
			}
			double value_next = -temp / mat.get_point(n, n);
			next.set_value_linear(n, value_next);
			// check if grid point changes more than current largest change
			if(abs(value_next - prev.get_value_linear(n)) > largest_change){
				largest_change = abs(value_next - prev.get_value_linear(n));
			}
		}
	}
	// fundamental breaking condition - effective equilibrium occurs when largest change is less than
	// the iteration tolerance level set by the user
	if(largest_change < tolerate){
		keep_iterating = false;
	}
	// vectors have been updated by reference, now return signal to timeloop
	return keep_iterating;
}


void timeloop(Stencil& stencil, MMatrix& matrix) {
	int maximum_time_in_seconds;
	float dt;
	// time is an integer, time step-size (dt) is a float
	// setting "dt" to 1 essentially makes "time" into "maximum iterations"
	cout << "Maximum time in seconds, and time step-size (dt) :" << endl;
	cin >> maximum_time_in_seconds >> dt >> transform;  
	cout << "Input your tolerance value" << endl;
	double tolerance; // iteration tolerance, for use in defining effective equilibrium
	cin >> tolerance;
	// make two linearized vectors
	Linear C(stencil);
	Linear D = C;
	bool C_next = true; // boolean for book-keeping
	bool keep_iterating = true; // booleance for cut-off when at effective equilibrium
	double final_time = 0.0; // tracking the final "time" in seconds
	for (double time = 0.0; time < maximum_time_in_seconds; time += dt) {
		C_next = !C_next;
		if (C_next) {
			keep_iterating = next_Jacobi(D, C, matrix, tolerance);
		}
		else {
			keep_iterating = next_Jacobi(C, D, matrix, tolerance);
		}
		final_time = time;
		if(!keep_iterating){
			break;
		}
	}
}
int main() {
	int n; // storing the problem number
	cout << "Which situation do you want to see? 1, 2, 3, or 4?" << endl;
	cin >> n; 
	Stencil stencil(n); // construct initial conditions object
	MMatrix matrix(stencil.get_values()); // construct "M" matrix
	timeloop(stencil, matrix); // start the timeloop
	return 0;
}