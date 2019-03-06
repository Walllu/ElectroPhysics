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
	Grid() {
	}
	Grid(int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg) : x_dim(x_dim_arg), y_dim(y_dim_arg), x_phys(x_phys_arg), y_phys(y_phys_arg) {
		// this is a more flexible default starter grid of zeros
		// just define the x,y array size and physical depth of each
		// and we set dx, dy here based on physical_length/grid_spaces
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
		return x_phys;
	}
	float get_y_phys() {
		return y_phys;
	}
	float get_x_position(int i) {
		return i * dx - x_phys / 2.0;
	}
	float get_y_position(int i) {
		return i * dy - y_phys / 2.0;
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

class MMatrix {
private:
	vector < vector<int> > M; //linearized potential vector
	int M_dim; // length of the arrays
	int x_dim; // x dim
	int y_dim; // y dim

public:
	MMatrix(Grid g) {
		x_dim = g.get_x_dim();
		y_dim = g.get_y_dim();
		M_dim = x_dim * y_dim;

		vector < vector<int> > temp(M_dim);

		for (int i = 0; i < x_dim; i++) {
			for (int j = 0; j < x_dim; j++) {
				int current_point = i + j * x_dim;
				vector<int> position;
				if (current_point < x_dim) {
					// current point is in first row
					position.push_back(current_point + x_dim);
					//Mrow[current_point + x_dim] = 1; // below
					//Mrow[M_dim - x_dim + current_point] = 1; // --------------------------------- pacman - first row
					position.push_back(current_point +M_dim - x_dim);
					if (current_point == 0) { // ----------------first row, first col
						//Mrow[current_point + 1] = 1;
						position.push_back(current_point + 1);
						//Mrow[x_dim - 1] = 1; //------------------------------------------ pacman - first row
						position.push_back(x_dim - 1);
					}
					else if (current_point == x_dim - 1) { //----first row, last col
						//Mrow[current_point - 1] = 1;
						position.push_back(current_point - 1);
						
						//Mrow[0] = 1; //------------------------------------------------ pacman - first row
						position.push_back(0);
					}
					else {
						position.push_back(current_point - 1);
						position.push_back(current_point + 1);
						//Mrow[current_point + 1] = 1;
						//Mrow[current_point - 1] = 1;
					}
				}
				else if (current_point >= M_dim - x_dim) {
					// current point is in last row
					position.push_back(current_point - x_dim);
					position.push_back(current_point + x_dim - M_dim);
					//Mrow[current_point - x_dim] = 1; // above
					//Mrow[current_point - M_dim + x_dim] = 1;//------------------------------- pacman - last row
					if (current_point == M_dim - x_dim) { //-------last row, first col
						position.push_back(current_point + 1);
						position.push_back(M_dim - 1);
						//Mrow[current_point + 1] = 1;
						//Mrow[M_dim - 1] = 1; //------------------------------------------------ pacman - last row
						//Mrow[0] = 1;//------------------------------------------------------ pacman - last row
					}
					else if (current_point == M_dim - 1) {//-------last row, last col
						position.push_back(current_point - 1);
						position.push_back(M_dim - x_dim);
						//Mrow[current_point - 1] = 1;
						//Mrow[M_dim - x_dim] = 1; //----------------------------------------- pacman - last row
						//Mrow[x_dim-1] = 1;//------------------------------------------------ pacman - last row
					}
					else {
						position.push_back(current_point + 1);
						position.push_back(current_point - 1);
						//Mrow[current_point + 1] = 1;
						//Mrow[current_point - 1] = 1;
					}
				}
				else if (current_point % x_dim == 0 || current_point % x_dim == x_dim - 1) {
					position.push_back(current_point + x_dim);
					position.push_back(current_point - x_dim);
					//Mrow[current_point + x_dim] = 1; // below
					//Mrow[current_point - x_dim] = 1; // above
					if (current_point % x_dim == 0) {
						// do something if in first column in last row
						position.push_back(current_point + 1);
						position.push_back(current_point - 1);
						//Mrow[current_point + 1] = 1;
						//Mrow[current_point + x_dim - 1] = 1; //----------------------------------- pacman - first col
					}
					else {
						position.push_back(current_point - 1);
						position.push_back(current_point + 1);
						//Mrow[current_point - 1] = 1;
						//Mrow[current_point - x_dim + 1] = 1; //----------------------------------- pacman - last col
					}

				}
				else {
					position.push_back(current_point + x_dim);
					position.push_back(current_point - x_dim);
					position.push_back(current_point + 1);
					position.push_back(current_point - 1);
					//Mrow[current_point + x_dim] = 1; // below
					//Mrow[current_point - x_dim] = 1; // above
					//Mrow[current_point + 1] = 1; // right
					//Mrow[current_point - 1] = 1; // left
				}

				
				temp[current_point] = position;
			}
		}

		M = temp;
	};
	float get_size() {
		return M_dim;
	}
	int get_x_dim() {
		return x_dim;
	}
	int get_y_dim() {
		return y_dim;
	}
	vector <int> get_list(int i) {
		return M[i];
	}
	
};


// class for making the initial value stencils for problems
class Stencil {
private:
	// Grids that represnet the locations and values of all the constants
	Grid problem_values, problem_constant_locations;
	// problem_values are the 'locations and values of all points in the initial moment
	// problem_constant_locations has the locations of all constant values in the grid
public:
	Stencil(int num, int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg) { // we make a constructor that acts differently depending on which input you start with
		//cout << "Input number of x, y dimensions and x, y distances" << endl;
		int x_dim = x_dim_arg;
		int y_dim = y_dim_arg;
		float x_phys = x_phys_arg;
		float y_phys = y_phys_arg;
		//cin >> x_dim >> y_dim >> x_phys >> y_phys;
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
		}
		else if (num == 2) {

			//create problem 1, I'll just do hot walls)

			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					if (j == initial_values.get_y_dim() - 1) {

						initial_values.set_point(i, j, 20);
						constant_locations.set_point(i, j, 1);
					}
					else if (j == 0) {
						initial_values.set_point(i, j, -20);
						constant_locations.set_point(i, j, 1);
					}
					else {

						initial_values.set_point(i, j, 0);
						constant_locations.set_point(i, j, 0);
					}
				}
			}
		}
		else if (num == 3) {

			//create problem 3, two cylinders)
			float x, y;
			float out_r = 5;
			float inn_r = 1;
			//cout << "type inner and outer radius" << endl;
			//cin >> inn_r >> out_r;
			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				x = initial_values.get_x_position(i);
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					y = initial_values.get_y_position(j);
					if (x*x + y * y >= out_r * out_r) {

						initial_values.set_point(i, j, 10);
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

		}
		else if (num == 4) {
			//create problem 4, cylinder between walls)
			float x, y;
			float inn_r = 1;

			for (int i = 0; i < initial_values.get_x_dim(); i++) {
				x = initial_values.get_x_position(i);
				for (int j = 0; j < initial_values.get_y_dim(); j++) {
					y = initial_values.get_y_position(j);
					if (j == 0) {

						initial_values.set_point(i, j, 10);
						constant_locations.set_point(i, j, 1);
					}
					else if (j == initial_values.get_y_dim() - 1) {

						initial_values.set_point(i, j, -10);
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
		}

		else if (num == 5) {

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
			double a = 1;
			double b = 0.5; // x y size of the rectangles

			double V = 10;

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
		// end of problem 1


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
		vector<int> make_vec_const(size, 0);

		for (int i = 0; i < x_dim; i++) {
			for (int j = 0; j < x_dim; j++) {
				temp_vec_potential[i + j * x_dim] = s.get_value_at_point(i, j);
				if (s.is_it_constant(i, j)) {
					make_vec_const[i + j * x_dim] = 1;
				}
				else {
					make_vec_const[i + j * x_dim] = 0;
				}

			}
		}

		vec_potential = temp_vec_potential;

		vec_const = make_vec_const;
	};
	float get_size() {
		return size;
	}
	int get_x_dim() {
		return x_dim;
	}
	int get_y_dim() {
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

// function for obtain next Linear class object from previous one
// modified now to look at the tolerance of the iterations - returning false will signal to the timeloop to cut iteration short
bool next_GaussSeidel(Linear& next, MMatrix& mat, double tolerate) {
	double temp;
	int size = next.get_size();
	double largest_change = 0;

	vector<int> position;

	bool keep_iterating = true; // ---------- this boolean will inform the time loop to keep iterating, given tolerance
	double prev_value; // going to use this to save the value of the linear vector at the point, in order to compare difference
	//--------------- The biggest difference between the Jacobi and the Gauss-Seidel methods is that the former only relies on point from the previous
	//--------------- iteration in time to construct the next one, while the latter also uses new points that have already been evaluated
	for (int n = 0; n < size; n++) {
		// if the current point is not constant...
		position = mat.get_list(n);
		if (!(next.is_it_constant_linear(n))) {
			temp = 0;
			prev_value = next.get_value_linear(n);
			for (int index = 0; index < 4; index++) {
				temp += next.get_value_linear(position[index]);
			}
			double value_next = -temp / (-4.0);
			next.set_value_linear(n, value_next);
			if (abs(value_next - prev_value) > largest_change) {
				largest_change = abs(value_next - prev_value);
			}
		}
	}
	if (largest_change < tolerate) {
		keep_iterating = false;
	}
	return keep_iterating;
}

vector< vector<double> > timeloop(Stencil& stencil, MMatrix& mat) {
	// working with default initial conditions dx, dy, dt, tmax
	int maximum_time_in_seconds = 1000; //transform;
	float dt = 0.01;
	// time is an integer, time step-size (dt) is a float, and transform is an integer used for debugging purposes - enter 0 for transform if you don'e give a shit about debugging
	//cout << "Maximum time in seconds, and time step-size (dt) : transform too" << endl;
	//cin >> maximum_time_in_seconds >> dt >> transform;  // "transform" is the iteration number that you want to view - transform = 1 is the first time step, 2 is 2nd, etc.
	//cout << "Input your tolerance value" << endl;
	double tolerance = 0.0001;
	//cin >> tolerance;
	// make one linearized vectors
	Linear C(stencil);
	// test print one
	// iterating through time using the linearized Jacobi method
	bool keep_iterating = true;
	double final_time = 0.0;
	for (double time = 0.0; time < maximum_time_in_seconds; time += dt) {
		keep_iterating = next_GaussSeidel(C, mat, tolerance);
		final_time = time;
		// now we check if the tolerance has been reached
		if (!keep_iterating) {
			break;
		}
	}

	int x_dim = C.get_x_dim();
	int y_dim = C.get_y_dim();
	vector< vector<double> > print;

	vector<double> v;
	for (int i = 0; i < x_dim; i++) {
		for (int j = 0; j < y_dim; j++) {
			v.push_back(C.get_value_linear(i + j * x_dim));
		}

		print.push_back(v);
		v.clear();
		//cout << "hey" << endl;
	}

	return print;
}



int main() {

	return 0;
}
