#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


// a helper class that defines the grid and grid spacing
class Grid {
protected:
	int x_dim, y_dim; // number of nodes in x, y direction, (x,y index goes from 0 to x_dim -1 , y_dim - 1)
	float x_phys, y_phys; // physical x, y distances
	float dx, dy;   // grid spacing in each x,y direction
	vector<vector <double> > grid; // the actual grid of some dimension
public:
	Grid() {
	}
	Grid(int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg) : x_dim(x_dim_arg), y_dim(y_dim_arg), x_phys(x_phys_arg), y_phys(y_phys_arg) {
		// this is the default starter grid of zeros
		// just define the x,y array size and physical length of each
		// and we set dx, dy here based on physical_length/grid_spaces
		dx = x_phys / (x_dim - 1);
		dy = y_phys / (y_dim - 1);
		// construct the grid of zeros
		vector <double> y_axis(y_dim);
		vector < vector<double> > x_axis;
		for (int i = 0; i < x_dim; i++) {
			x_axis.push_back(y_axis);
		}
		grid = x_axis;
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
		return x_phys; // getter for physical length in x dimension
	}
	float get_y_phys() {
		return y_phys; // getter for physical length in y dimension
	}
	float get_x_position(int i) {
		return i * dx - x_phys / 2.0; // // calculate x position with respect to central point depending on index
	}
	float get_y_position(int i) {
		return i * dy - y_phys / 2.0; // calculate y position with respect to central point depending on index
	}
	void set_point(int i, int j, double val) {
		grid[i][j] = val;	// sets a point in the grid to "val"
							// i refers to the x position, j to the y position
	}
	double get_point(int i, int j) {
		return grid[i][j];	// getter for value of point in grid
							// i refers to the y position, j to the x position
	}
};

// helper class representing the matrix
class MMatrix {
private:
	vector < vector<int> > M; //linearized potential vector
	int M_dim; // dimension of the matrix
	int x_dim; // number of nodes in x direction
	int y_dim; // number of nodes in y dim

public:
	MMatrix(Grid g) {
		x_dim = g.get_x_dim(); // get number of nodes in x direction
		y_dim = g.get_y_dim(); // get number of nodes in y direction
		M_dim = x_dim * y_dim; // calculate the dimension of the matrix

		//vector with entries corresponding to the positions on the potential vector
		//at each position there is a vector of potential's positions that need to be added to calculate the nabla operator
		vector < vector<int> > temp(M_dim);

		// for each point on the grid create vector with positions that need to be added
		for (int i = 0; i < x_dim; i++) {
			for (int j = 0; j < x_dim; j++) {
				int current_point = i + j * x_dim; // ---------------------------------------------linearize grid indexes
				vector<int> position; // --------------------------------------------temporary vector to write indexes to
				if (current_point < x_dim) { //------------------------- first row
					position.push_back(current_point + x_dim); // ------ below 
					position.push_back(current_point + M_dim - x_dim); //above pacman to last row
					if (current_point == 0) { // ------------------------first row, first col
						position.push_back(current_point + 1); // ------ right
						position.push_back(x_dim - 1); //--------------- left pacman to last column
					}
					else if (current_point == x_dim - 1) { //------------ first row, last col
						position.push_back(current_point - 1); // ------- left
						position.push_back(0); // ----------------------- right pacman to first column
					}
					else {
						position.push_back(current_point - 1); // ------- left
						position.push_back(current_point + 1); // ------- right
					}
				}
				else if (current_point >= M_dim - x_dim) { // ----------- last row
					position.push_back(current_point - x_dim); // ------- above
					position.push_back(current_point + x_dim - M_dim); // below pacman to first row
					if (current_point == M_dim - x_dim) { //------------- last row, first col
						position.push_back(current_point + 1); // ------- right
						position.push_back(M_dim - 1); // --------------- left pacman to first column
					}
					else if (current_point == M_dim - 1) {//------------- last row, last col
						position.push_back(current_point - 1); // ------- left
						position.push_back(M_dim - x_dim); // ----------- right pacman to last column
					}
					else {
						position.push_back(current_point + 1); // ------- right
						position.push_back(current_point - 1); // ------- left
					}
				}
				else if (current_point % x_dim == 0 || current_point % x_dim == x_dim - 1) {
					// -------------------------------------------------- first or last column (not first or last row)
					position.push_back(current_point + x_dim); // ------- below
					position.push_back(current_point - x_dim); // ------- above
					if (current_point % x_dim == 0) { // ---------------- first column
						position.push_back(current_point + 1); // ------- right
						position.push_back(current_point + x_dim - 1); // left pacman to last column
					}
					else { // ------------------------------------------- last column
						position.push_back(current_point - 1); // ------- left
						position.push_back(current_point - x_dim + 1); // right pacman to first column
					}
				else { // --------------------------------------- --------point not at the edge of the grid
					position.push_back(current_point + x_dim); // ------- below
					position.push_back(current_point - x_dim); // ------- above
					position.push_back(current_point + 1); // ----------- right
					position.push_back(current_point - 1); // ----------- left
				}


				temp[current_point] = position; // assign list of positions to potential at vector
				}
			}
			M = temp; // set matrix vector
		};
		float get_size() {
			return M_dim; // getter for matrix size
		}
		int get_x_dim() {
			return x_dim; // getter for grid x dimension
		}
		int get_y_dim() {
			return y_dim; // getter for grid y dimension
		}
		vector <int> get_list(int i) { // get list of positions for particular potential point
			return M[i];
		}

	};

	// class for making the initial condition for problems
	class Stencil {
	private:
		// Grids that represnet the locations and values of all the constants
		Grid problem_values, problem_constant_locations;
		// problem_values are the 'locations and values of all points in the initial moment
		// problem_constant_locations has the locations of all constant values in the grid
	public:
		// we make a constructor that acts differently depending on which input you start with
		Stencil(int num, int x_dim_arg, int y_dim_arg, float x_phys_arg, float y_phys_arg) {
			int x_dim = x_dim_arg;
			int y_dim = y_dim_arg;
			float x_phys = x_phys_arg;
			float y_phys = y_phys_arg;
			// create grids temporary grids for potential values and constant locations
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

				//create problem 2, we will just do one wall hot and other cold)

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
				float x, y; // x and y physical positions
				float out_r = 5; // outer radius
				float inn_r = 1; // inner radius
				for (int i = 0; i < initial_values.get_x_dim(); i++) {
					x = initial_values.get_x_position(i);
					for (int j = 0; j < initial_values.get_y_dim(); j++) {
						y = initial_values.get_y_position(j);
						if (x*x + y * y >= out_r * out_r) { // check if not within large cylinder

							initial_values.set_point(i, j, 10);
							constant_locations.set_point(i, j, 1);
						}
						else if (x*x + y * y <= inn_r * inn_r) { // check if within small cylinder

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 1);
						}
						else { // between cylinders

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 0);
						}
					}
				}

			}
			else if (num == 4) {
				//create problem 4, cylinder between plates)
				float x, y; // x and y physical positions
				float inn_r = 1; // radius

				for (int i = 0; i < initial_values.get_x_dim(); i++) {
					x = initial_values.get_x_position(i);
					for (int j = 0; j < initial_values.get_y_dim(); j++) {
						y = initial_values.get_y_position(j);
						if (j == 0) { // check if at the plate position

							initial_values.set_point(i, j, 10);
							constant_locations.set_point(i, j, 1);
						}
						else if (j == initial_values.get_y_dim() - 1) { // check if at the plate position

							initial_values.set_point(i, j, -10);
							constant_locations.set_point(i, j, 1);
						}
						else if (x*x + y * y <= inn_r * inn_r) { // check if within cylinder

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 1);
						}
						else { // between plates and cylinder

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 0);
						}
					}
				}
			}

			else if (num == 5) {

				//create problem 5, edge coupled stripline

				for (int i = 0; i < initial_values.get_x_dim(); i++) {
					for (int j = 0; j < initial_values.get_y_dim(); j++) {
						if (j == initial_values.get_y_dim() - 1) { // check if at plate position

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 1);
						}
						else if (j == 0) { // check if at plate position
							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 1);
						}
						else { // rest of points

							initial_values.set_point(i, j, 0);
							constant_locations.set_point(i, j, 0);
						}
					}
				}

				// x y size of the rectangles
				double a = 1;
				double b = 0.5;

				double V = 10; // absolute value of potentials of the conductors

				for (int i = 0; i < initial_values.get_x_dim(); i++) {
					for (int j = 0; j < initial_values.get_y_dim(); j++) {
						if (abs(initial_values.get_y_position(j)) <= b / 2.0) {
							// check if at positive constant conductor
							if (abs(initial_values.get_x_position(i) - initial_values.get_x_phys() / 8.0) <= a / 2.0) {
								initial_values.set_point(i, j, V);
								constant_locations.set_point(i, j, 1);
							}
							// check if at negative constant conductor
							else if (abs(initial_values.get_x_position(i) + initial_values.get_x_phys() / 8.0) <= a / 2.0) {
								initial_values.set_point(i, j, -V);
								constant_locations.set_point(i, j, 1);
							}
							// check if at first grounded constant conductor
							else if (abs(initial_values.get_x_position(i) - 3.0*initial_values.get_x_phys() / 8.0) <= a / 2.0) {
								initial_values.set_point(i, j, 0);
								constant_locations.set_point(i, j, 1);
							}
							// check if at second grounded constant conductor
							else if (abs(initial_values.get_x_position(i) + 3.0*initial_values.get_x_phys() / 8.0) <= a / 2.0) {
								initial_values.set_point(i, j, 0);
								constant_locations.set_point(i, j, 1);
							}
						}
					}
				}

			}

			problem_values = initial_values;
			problem_constant_locations = constant_locations;


		}

		Grid get_values() {
			return problem_values; // getter for grid of potential values
		}

		Grid get_constants() {
			return problem_constant_locations; // getter for grid indicationg constants
		}

		double get_value_at_point(int i, int j) {
			return problem_values.get_point(i, j); // getter for potential value at point of grid
		}

		bool is_it_constant(int i, int j) { // check if ptoential is constant at point on the grid
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
			x_dim = s.get_values().get_x_dim(); // get x grid dimension
			y_dim = s.get_values().get_y_dim(); // get y grid dimension
			size = x_dim * y_dim; // calculate length of linearized potential vector

			vector<double> temp_vec_potential(size, 0); // vector to store the potential values
			vector<int> make_vec_const(size, 0); // create vector to indicate constants

			// assign potentials and constants form grid to vector
			for (int i = 0; i < x_dim; i++) {
				for (int j = 0; j < x_dim; j++) {
					temp_vec_potential[i + j * x_dim] = s.get_value_at_point(i, j); // index in vector is i + x_dim * j
					if (s.is_it_constant(i, j)) {
						make_vec_const[i + j * x_dim] = 1;
					}
					else {
						make_vec_const[i + j * x_dim] = 0;
					}

				}
			}

			vec_potential = temp_vec_potential; // set vector storing potential

			vec_const = make_vec_const; // set vector storing constants
		};
		float get_size() {
			return size; // getter for length of vector
		}
		int get_x_dim() {
			return x_dim; // getter for x drid dimension
		}
		int get_y_dim() {
			return y_dim; // getter for y grid dimension
		}
		double get_value_linear(int i) {
			return vec_potential[i]; // return potential value at position on linerized vector 
		}
		void set_value_linear(int i, double val) {
			vec_potential[i] = val; // set potential value at position on linerized vector
		}


		bool is_it_constant_linear(int i) { // check if potential is constant at the position
			if (vec_const[i] == 1) {
				return true;
			}
			else {
				return false;
			}
		}
	};

	// function iterativery applied to potential vector to calculate the next values
	// retyrns true if precision is not reached and false if it is reached
	bool next_SOR(Linear& next, MMatrix& mat, double tolerate, double sor) {
		double temp; // temporarily sums values of potential at the locations
		int size = next.get_size(); // length of potential vector
		double largest_change = 0; // largest change in single iteration watch

		vector<int> position; // vector to read list of indexes from MMatrix class for each point

		bool keep_iterating = true; // ---------- this boolean will inform the time loop to keep iterating, given tolerance
		double prev_value; // going to use this to save the value of the linear vector at the point, in order to compare difference
		//--------------- The biggest difference between the Jacobi and the Gauss-Seidel methods is that the former only relies on point from the previous
		//--------------- iteration in time to construct the next one, while the latter also uses new points that have already been evaluated
		for (int n = 0; n < size; n++) { // loop the vector of potential 
			// if the current point is not constant...
			position = mat.get_list(n); // read list of indexes
			if (!(next.is_it_constant_linear(n))) { // check if potential at point is not constant
				temp = 0;
				prev_value = next.get_value_linear(n); // store potential value at the point
				for (int index = 0; index < 4; index++) { // add potentials at four positions
					temp += next.get_value_linear(position[index]);
				}
				double value_next = prev_value * (1 - sor) - sor * temp / (-4.0); // calculate next value of potential
				next.set_value_linear(n, value_next); // update potential value
				if (abs(value_next - prev_value) > largest_change) { // update the largest change
					largest_change = abs(value_next - prev_value);
				}
			}
		}

		if (largest_change < tolerate) { //------------------------------------ check if largest change is below our toleration
			keep_iterating = false;
		}
		return keep_iterating;
	}


	vector< vector<double> > timeloop(Stencil& stencil, MMatrix& matrix) {
		// returns final grid of potential values
		int maximum_time_in_seconds = 1000;
		float dt = 0.01;
		// time is an integer, time step-size (dt) is a float
		// serve to count the iteratins done

		double tolerance = 0.0001; // set tolerance

		// make one linearized vector of potential with initial values
		Linear C(stencil);
		// calculate successive overrelaxation parameter
		double sor = 2. / (1 + sin(3.14159265358979323846 / (stencil.get_values().get_x_dim())));

		// iterating through time using the SOR method
		bool keep_iterating = true;
		double final_time = 0.0;
		for (double time = 0.0; time < maximum_time_in_seconds; time += dt) {
			keep_iterating = next_SOR(C, matrix, tolerance, sor);
			final_time = time;
			// now we check if the tolerance has been reached
			if (!keep_iterating) {
				break;
			}
		}

		// transform linear potential into grid
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
		}

		return print;
	}



	int main() {

		return 0;
	}