# ode_solver: about
A dynamical library for solving the systems of Ordinary Differential Equations (ODEs) of the first order. Three basic schemes are implemented: Euler scheme; Runge-Kutta scheme of the 4th order; Gear scheme (a.k.a. BDF method) of the 2nd, 4th and 6th order.

# Equations
The library can be used to solved the system of the ODEs:

$$
\begin{align}
  \frac{d\vec{y}}{dt} &=& \vec{f}(t,\,\vec{y}),\\
  \vec{y}(t=0) &=& \vec{y}_0,
\end{align}
$$

where $\vec{y} = [y^{(0)}, y^{(1)}, y^{(2)}, ... y^{(N-1)}]^T$ - vector (column) of unknown variables, $N$ - number of equations (variables)  $\vec{f}(t, \vec{y}) = [f^{(0)}, f^{(1)}, f^{(2)}, ... f^{(N-1)}]^T$ - vector of the derivatives in the system (a.k.a _right-hand side_ of the system), $\vec{y}_0 = [y^{(0)}_0, y^{(1)}_0, y^{(2)}_0, ... y^{(N-1)}_0]^T$ - corresponding vector of the initial values.

# Code
## Structure
The code is organised into the following directories:
- `./lib/`: source files of the library
- `./test_src/`: source files for testing the library on the problems with known exact solution
- `./cb_project/`: Code::Blocks project files to edit and compile the code
- `./msvs_project/`: Microsoft Visual Studio project files to edit and compile the code

The library source code consists of the following files:
- `utypes.h/.cpp`: user defined types and function to handle them
- `ode.h/.cpp`: definition and implementation of the ODE class used to represent the equations to be solved
- `ode_solver.h/.cpp`: definition and implementation of the Basic ODE Solver class. The class implements common functions for any ODE solver and used as a basic class to create corresponding classes for specific schemes
- `euler_solver.h/.cpp`: definition and implementation of the Euler ODE Solver based on the basic ODE Solver Class
- `rk_solver.h/.cpp`: definition and implementation of the Runge-Kutta ODE Solver based on the basic ODE Solver Class
- `gear_solver.h/.cpp`: definition and implementation of the Gear ODE Solver based on the basic ODE Solver Class

The source code for testing the library:
- `main.cpp`: entry point for the program testing the library schemes
- `ctest.h/.cpp`: definition and implementation of the basic Test Class. The class allows to set specific equations to be solved, pick a a scheme to be used, run the solution process and save the results.
- `test1.h/.cpp`: definition and implementation of the Test 1 based on the basic Test Class.
- `test2.h/.cpp`: definition and implementation of the Test 1 based on the basic Test Class.
- `test3.h/.cpp`: definition and implementation of the Test 1 based on the basic Test Class.
- `test4.h/.cpp`: definition and implementation of the Test 1 based on the basic Test Class.
- `test5.h/.cpp`: definition and implementation of the Test 1 based on the basic Test Class.

## Class `TODE`
Methods description

## Class `TRungeKuttaSolver`
Methods description

# Usage
Below are the steps to configure the library for solving a specific problem and to make one step of integration.

1) Add the library to your project by including corresponding headers. E.g., for solution of the ODE via the Runge-Kutta scheme:
```
#include "./lib/ode.h"
#include "./lib/rk_solver.h"
```
2) Define model parameters
```
/// equation parameters

size_t eqs_number = 1; // number of equations
double y0 = 1.0;       // initial condition for a single equation
double tau = 1.0;      // typical time scale for a specific problem

/// solver parameters

/// tolerance - required relative accuracy in the scheme with step size control
double tol = 1e-4;
/// initial time step
double dt = 1e-4 * tau;
/// minimum allowed time step in the scheme with step size control
double dt_min = 1e-6 * tau;
/// maximum allowed time step in the scheme with step size control
double dt_max = 5e-1 * tau;
/// the order of the scheme
size_t order = 4;
/// flag: is solution found (true) or not (false)
bool is_solution_found = false;

/// Set and configure model parameters

double t = 0.0;          // current time
std::vector<double> ic;  // vector of initial conditions
std::vector<double> sol; // vector of the solution after one step of integration

// add the initial value to the initial values vector (example for a single equation)
ic.push_back(y0);

// allocate memory for the solution vector
sol.resize(eqs_number);
```

2) Create instances of the ODE class (type: `TODEs`) and solver (type: `TRungeKuttaSolver`/`TEulerSolver`/`TGearSolver`). E.g.,
```
/// object to handle the system of ODEs
TODEs ode;
/// object to handle the scheme (solver)
TRungeKuttaSolver solver;
```
3) Configure object 'ode' (equation/system of equation).
Let the RHS of the equation is defined as
```
double dydt(const double t, std::vector<double> y);
```

Exact solution:
```
double y_exact(const double t);
```

To add the system to the object 'ode' one need to define the pointers to the RHS function and exact solution, then import these pointers into the 'ode' object using the interface functions:
```
/// pointer to a RHS function
PFunctionVectArg rhs = dydt;
/// pointer to an exact solution
PFunction pSol = &y_exact;
/// add the equation to the object 'ode'
ode.AddRHS(rhs);
/// add corresponding initial condition to the object 'ode'
ode.AddIC(y0);
/// add the exact solution to the object 'ode'
ode.AddSolution(pSol);
```

Any additional equations of the system can be added in a similar way.

4) Configure the solver.
```
// set the equations number in the solver
solver.SetVariablesNumber(eqs_number);
// set the order of the scheme
solver.SetOrder(order);

// allocate internal variables in the solver
solver.AllocateSolution();
solver.AllocateCoefficients();

// activate step size control
solver.ActivateStepSizeControl();
// set required relative accuarcy (tolerance)
solver.SetTolerance(tol);
// set maximum allowed time step
solver.SetMaxStep(dt_max);
// set minimum allowed time step
solver.SetMinStep(dt_min);
// import the system of ODEs into solver
solver.ImportODEs(&ode);
```

5) Configure initial state of the problem:
```
// set inital time to zero
solver.SetTime(0);
// set initial time step
solver.SetTimeStep(dt);
// import the initial solution into the solver
solver.ReadInitialVector(sol, eqs_number);
```

6) Run the integration and get the solution
```
// make an integration step
is_solution_found = solver.Solve();

if (!is_solution_found)
  printf("Warning! Solution was not found. Exiting.\n");
else
{
  // write the solution into the solution vector

  for (size_t i_sol = 0; i_sol < eqs_number; i_sol++)
    sol.at(i_sol) = solver.GetSolution(i_sol);

  // prepare solver to a next time step

  dt = solver.GetTimeStep(); // get the update time step
  solver.UpdateStep(dt);     // write the updated time step into the solver
  t += dt;                   // update time
  solver.SetTime(t);         // set the updated time in the solver
}
```

## Example code
Below is the example code for solving a single ODE using the RK scheme of the 4th order with step size control.
```
#include <iostream>
#include <cmath>
#include "./lib/ode.h"
#include "./lib/rk_solver.h"

using namespace std;

/// function defining the derivative in the ODE
double dydt(const double t, std::vector<double> y);
/// exact solution of the ODE
double y_exact(const double t);

int main()
{
    cout << "Hello! This program tests the ODE solver" << endl;
    cout << endl;

    /// equation parameters

    size_t eqs_number = 1; // number of equations
    double y0 = 1.0;       // initial condition for a single equation
    double tau = 1.0;      // typical time scale for a specific problem

    /// solver parameters

    /// tolerance - required relative accuracy in the scheme with step size control
    double tol = 1e-4;
    /// initial time step
    double dt = 1e-2 * tau;
    /// minimum allowed time step in the scheme with step size control
    double dt_min = 1e-6 * tau;
    /// maximum allowed time step in the scheme with step size control
    double dt_max = 5e-1 * tau;
    /// the order of the scheme
    size_t order = 4;
    /// flag: is solution found (true) or not (false)
    bool is_solution_found = false;

    /// Set and configure model parameters

    double t = 0.0;          // current time
    std::vector<double> ic;  // vector of initial conditions
    std::vector<double> sol; // vector of the solution after one step of integration

    ic.push_back(y0);

    // allocate memory for the solution vector
    sol.resize(eqs_number);

    /// 2) create ODE and solver instances

    TODEs ode;
    TRungeKuttaSolver solver;

    /// 3) Configure ODE

    /// pointer to a RHS function
    PFunctionVectArg rhs = dydt;
    /// pointer to an exact solution
    PFunction pSol = &y_exact;
    /// add the equation to the object 'ode'
    ode.AddRHS(rhs);
    /// add corresponding initial condition to the object 'ode'
    ode.AddIC(y0);
    /// add the exact solution to the object 'ode'
    ode.AddSolution(pSol);

    /// 4) Configure solver

    // set the equations number in the solver
    solver.SetVariablesNumber(eqs_number);
    // set the order of the scheme
    solver.SetOrder(order);


    // allocate internal variables in the solver
    solver.AllocateSolution();
    solver.AllocateCoefficients();

    // activate step size control
    solver.ActivateStepSizeControl();
    // set required relative accuarcy (tolerance)
    solver.SetTolerance(tol);
    // set maximum allowed time step
    solver.SetMaxStep(dt_max);
    // set minimum allowed time step
    solver.SetMinStep(dt_min);
    // import the system of ODEs into solver
    solver.ImportODEs(&ode);

    /// 5) Configure initial state of the problem

    // set inital time to zero
    solver.SetTime(0);
    // set initial time step
    solver.SetTimeStep(dt);
    // import the initial solution into the solver
    solver.ReadInitialVector(ic, eqs_number);

    /// 6) Run the integration and get the solution

    // make an integration step
    is_solution_found = solver.Solve();

    if (!is_solution_found)
        printf("Warning! Solution was not found. Exiting.\n");
    else
    {
        // write the solution into the solution vector

        for (size_t i_sol = 0; i_sol < eqs_number; i_sol++)
            sol.at(i_sol) = solver.GetSolution(i_sol);

        // prepare solver to a next time step

        dt = solver.GetTimeStep(); // get the update time step
        solver.UpdateStep(dt);     // write the updated time step into the solver
        t += dt;                   // update time
        solver.SetTime(t);         // set the updated time in the solver

        // print numerical and exact solutions to standard output

        cout.setf(std::ios::scientific);
        cout.width(10);
        cout.precision(5);
        cout << "solution at t   = " << t << " is:" << endl;
        cout << "            y   = " << sol.at(0) << endl;
        cout << "exact solution:" << endl;
        cout << "            y_e = " << y_exact(t) << endl;
    }

    return 0;
}

double dydt(const double t, std::vector<double> y)
{
    return -10*y[0];
};

double y_exact(const double t)
{
    return exp(-10*t);
};
```

# Tests
