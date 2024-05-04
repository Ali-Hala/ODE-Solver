import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy import symbols, lambdify, simplify
import algebra

# This is the 2NSpace file. This is the main file that solves ODEs in the 2N space. In order to better understand this code, start reading at __main__ found on line 150.

# Define the ODE system with user-specified number of variables
def model(t, y, ode_exprs):
    dydt = [ode_expr(y, t) for ode_expr in ode_exprs]
    return dydt

def get_coefficient(row, indices, dimensions, iterations = 0):
    iterations += 1
    if iterations < dimensions:
        return get_coefficient(row[indices[iterations-1]-1], indices, dimensions, iterations)
    else:
        return row[indices[iterations - 1]-1]

def make_expr(index, init_row, iterations, dim, indices, terms, y_symbols, Original_K):
    iterations += 1 # Every time make_expr is recursively called, we jump into a higher dimension of the K tensor. The iterations variable is used to keep track of how many times we have jumped.
    for i, val in enumerate(index): # The index represents the element from K that we are working with.
        indices.append(i) # indices is essentially a list that will hold the address to a particular element in K. Every time we append a number to it, we add a coordinate to the address. We need a coordinate for each dimension of K.
        if iterations < dim: # If the number of dimensions we have jumped to is less than the total number of dimensions in K.
            make_expr(val, init_row, iterations, dim, indices, terms, y_symbols, Original_K) # Recursively call make_expr(), This time with an element of index being the argument for index.
        else: # else if iterations is equal to the total number of dimensions in the K tensor (i.e. we have the correct number of coordinates in indices).
            total_term = "" # total_term will be a string containing the y values that will be multiplied in a term.
            total_indices = [init_row] + indices # init_row is one more coordinate that needs to be added to the address of our element in K. It is essentially the coordinate for the first dimension.
            for _index in indices:
                total_term += str(y_symbols[_index]) + "*" # form a string of y_symbols multiplied with each other.
            total_term = total_term[:-1] # Remove the extra star at the end of the string.
            terms.append(str(get_coefficient(Original_K, [j + 1 for j in total_indices], dim, 0)) + "*" + total_term) # Get the number we want from K using the address in indices and the get_coefficient() function defined on line 14, then multiply it with total_term and append to terms.
        indices.pop() # Remove the last number in indices when we get to the end of the loop.

def get_constants():
    constants = {}
    while True:
        constant_input = input("Enter a constant (name=value), or press Enter to finish. You can only use Letters in constant names except for t, y, x, P, and d: ")
        if not constant_input:
            break
        try:
            name, value = constant_input.split("=")
            constants[name.strip()] = float(value.strip())
        except ValueError:
            print("Invalid constant format. Please use name=value.")
    return constants

def change_constants(constants):
    while True:
        print("Current Constants:")
        for name, value in constants.items():
            print(f"{name} = {value}") # Show the constants they currently exist.
        constant_name = input("Enter the name of the constant you want to change, or press Enter to finish: ")
        if not constant_name:
            break
        if constant_name in constants:
            new_value = input(f"Enter the new value for {constant_name}: ")
            try:
                constants[constant_name] = float(new_value)
                print(f"{constant_name} updated.")
            except ValueError:
                print("Invalid value. Please enter a number.")
        else:
            print("Constant not found.")

def replace_constants(equation, constants):
    for constant, value in constants.items():
        equation = equation.replace(constant, str(value))
    return equation

def solve_ode(equation_str, num_solutions, initial_guess, increments, equation_str2, resolution, num_variables, equation_str3):

    # Call the solve_equation() function from the algebra file. This function returns our Galerkin matrix (K) and its dimension.
    K, dimension = algebra.solve_equation(equation_str, num_solutions, initial_guess, increments, equation_str2, resolution, num_variables, equation_str3)

    # Getting the number of elements in the first dimension of the galerkin matrix. This number is the same as num_solutions on line 159
    first_layer = 0
    for i in K:
        first_layer += 1


    print("Dimension: " + str(dimension))
    print("Coefficients (K):", K)

    # Define our list of lambdified ODE expressions, our symbols for those expressions, and a list of those expressions as strings.
    ode_exprs = []
    y_symbols = symbols([f'y{j + 1}' for j in range(2 * first_layer)])
    t = symbols('t')
    expr_strs = []

    for i in range(first_layer): # Since we are working in the 2N space, this loop defines and lambdifies the first N ODEs. It then appends each ODE to the ode_exprs list.
        expr_strs.append(f'y{first_layer + i + 1}')
        ode_expr = lambdify((y_symbols, t), f'y{first_layer + i + 1}', 'numpy')
        ode_exprs.append(ode_expr)

    for i, val in enumerate(K): # This loop defines the rest of the ODEs in our system using the galerkin matrix.
        terms = [] # This list holds the terms for each ODE. It will reset to an empty list each time we define a new ODE (i.e. each time we go through the loop).
        iterations = 1 # The iterations variable is what is used to ensure that we index the correct number of dimensions when performing recursion of the make_expr() function. This may make more sense when looking at the definition of the make_expr() function on line 21.
        if dimension > 1:
            make_expr(val, i, iterations, dimension, [], terms, y_symbols, K)  # Calling the make_expr() function defined on line 21 if K has more than one dimension.
        else:
            terms.append(str(K[i]) + "*" + str(y_symbols[i])) # If K only has one dimension, then multiply each value in K with its corresponding y-value.
        str_together = "" # This string will contain the combination of all the terms in the ODE (i.e. the entire equation).
        for j in terms:
            str_together += str(j) + "+" # Combine all terms of the equation, each seperated by a '+' sign.
        expr_str = str_together[:-1] # Remove the extra '+' sign at the end.
        expr_strs.append(expr_str)
        ode_expr = lambdify((y_symbols, t), expr_str, 'numpy') # Lambdify the ODE.
        ode_exprs.append(ode_expr) # Add the ODE to the ode_exprs list.

    # Print generated ODE expressions
    print("ODE Expressions:")
    for i, expr in enumerate(expr_strs):
        print(f"y{i + 1}:", expr)

    # Set initial conditions
    y0 = [float(input(f"Enter initial condition for y{i + 1}: ")) for i in range(2 * first_layer)]

    end = float(input("Enter the end:")) # Set the end of the domain to be evaluated for and plotted.

    # Get user input for the resolution (number of time points)
    time_points = int(input("Enter the desired resolution (number of time points): "))

    # Solve the ODEs using solve_ivp with the desired resolution
    solution = solve_ivp(model, (0, end), y0, args=(ode_exprs,), method='RK45', t_eval=np.linspace(0, end, time_points),
                         atol=1e-8, rtol=1e-8)

    # Print the vectors corresponding to each variable and t
    print("Time values (t):", solution.t)

    for i, val in enumerate(solution.y):
        print(f"Solution values (y{i + 1}):", val)

    entry_count = 0
    for i in solution.t:
        entry_count += 1

    # Plot the results
    for i, val in enumerate(solution.y):
        plt.plot(solution.t, val, label=f'y{i + 1}')

    print("Number of entries: " + str(entry_count))
    plt.xlabel('Time')
    plt.ylabel('y(t)')
    plt.title('Numerical Integration of ODEs')
    plt.legend()
    plt.show()

if __name__ == "__main__":

    # The get_constants() function asks the user for any constants they want to use and returns a dictionary of those constants. The function is defined on line 36
    constants = get_constants()

    # This gets the user input for the characteristic equation. The equation is always 0 = (expression inputted by the user). Typically, this will be a trigonometric/hyperbolic-trigonometric equation.
    equation_str = input("Enter the expression for the characteristic equation in terms of y (use sin, cos, sinh, cosh): ")

    # Get user input for the number of solutions they want the computer to solve for the above equation.
    num_solutions = int(input("Enter the number of solutions you want: "))

    # Get user input for the initial guess since the computer will use the newton-raphson method to calculate solutions.
    initial_guess = float(input("Enter the initial guess (for Newton-Raphson): "))

    # This is the number of increments that the user wants to increment the initial guess by each time the computer calculates another solution.
    increments = float(input("Enter increments between each initial guess for each solution: "))

    # This is where the user inputs the expression for the mode shape in terms of x and y.
    equation_str2 = input("Enter right side of the mode shape in terms of x and y: ")

    # This is essentially the user telling the computer how many solutions for the mode shape they want.
    resolution = input("Enter desired resolution for the mode shape: ")


    # This is where the user inputs the number of dimensions they want in the Galerkin matrix.
    num_variables = int(input("Enter the number of desired dimensions for the Galerkin Matrix: "))

    # Here we ask the user to input the Galerkin integrand.
    equation_str3 = input(f"Enter an expression in terms of {', '.join([f'P{j+1}' for j in range(num_variables)])}: " )

    # After Everything is finished, we ask the user if they want to alter the values of any constants from the constants dictionary.
    running = True
    while running:
        again = int(input("Enter 1 to alter constant values or 0 to exit: "))
        if again == 0:
            running = False
        else:
            change_constants(constants) # Call the change_constants function defined on line 49. This function asks the user for the new values of the constants they have and changes them.

            # Solve everything again.
            solve_ode(replace_constants(equation_str, constants), num_solutions, initial_guess, increments, replace_constants(equation_str2, constants), resolution, num_variables, replace_constants(equation_str3, constants))
