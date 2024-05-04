# This is the file responsible for the creation of the Galerkin matrix used in the 2NSpace file
from scipy.optimize import fsolve
import numpy as np
from sympy import *
# Below is a class definition, followed by some function definitions. In order to better understand the code, start reading at main found on line 181.
class PofY: # This is the PofY class definition. The purpose of this class is to create objects that are similar to lists, only these lists behave a little diferently than normal lists.
    def __init__(self, vector, dvector=None, ddvector=None, dddvector=None, ddddvector=None): # This class has 5 attributes. These attributes are 5 lists.
        self.vector = vector # Please note that the term 'vector' in this case is simply used to refer to a list.
        self.dvector = dvector if dvector is not None else []
        self.ddvector = ddvector if ddvector is not None else []
        self.dddvector = dddvector if dddvector is not None else []
        self.ddddvector = ddddvector if ddddvector is not None else [] # When we create objects of this class on line 176, vector will be a row in matrix, dvector will be the corresponding row in dmatrix, and so on and so forth.
    # Lines 14-42 Deal with operator overloading. They are fairly simple since all they do is element-wise operation.
    def __mul__(self, other): # Element-wise multiplication.
        new = PofY([])
        for i, val in enumerate(self.vector):
            new.vector.append(val * other.vector[i])
        return new

    def __pow__(self, power, modulo=None): # Element-wise exponentiation.
        new = PofY([])
        for i in self.vector:
            new.vector.append(pow(i, power))
        return new

    def __add__(self, other): # Element-wise addition.
        new  = PofY([])
        for i, val in enumerate(self.vector):
            new.vector.append(val + other.vector[i])
        return new

    def __sub__(self, other): # Element-wise subtraction.
        new = PofY([])
        for i, val in enumerate(self.vector):
            new.vector.append(val - other.vector[i])
        return new

    def __truediv__(self, other): # Element-wise division.
        new = PofY([])
        for i, val in enumerate(self.vector):
            new.vector.append(val / other.vector[i])
        return new

    def integ(self, interval_width): # This function takes a PofY object and returns the integral of its vector using the trapezoidal method.
        result = 0
        for i,val in enumerate(self.vector):
            if i == 0 or i == len(self.vector)-1:
                result += val * (interval_width/2) # The interval_width is the width of each interval between each element in the vector/list. This length is defined on line 161.
            else:
                result += 2*val*(interval_width/2)
        return result

    def d(self, depth): # This function takes a PofY object and changes its vector to be one of its derivatives up to the fourth derivative.
        match depth:
            case 1:
                self.vector = self.dvector
            case 2:
                self.vector = self.ddvector
            case 3:
                self.vector = self.dddvector
            case 4:
                self.vector = self.ddddvector
        return self

def equation(x, expression):
    return lambdify(x, expression, "numpy")

def matrixify(matrix, equation, roots, x_span):
    for i in roots:
        inner_list = [] # The computer creates an inner list for each y value in the roots list.
        for j in x_span:
            inner_list.append(equation(j, i)) # The computer solves the equation for each x value in x_span using i as the y value and inserts each solution into the inner list.
        matrix.append(inner_list) # Add the inner list to the matrix, and repeat for each value in roots.
    print(matrix)

def store(array, indices, eq, dimension, iterations, last_row, interval_width):
    iterations += 1 # Increment iterations by 1 each time store() is recursively called, since iterations tells us how many dimensions we have jumped to.
    for i in array: # Iterate through array again, which in our case is matrix2.
        indices.append(i) # Append the current element in array to indices.
        if iterations < dimension:
            new_row = [] # If iterations is still < dimension, we want to put a new list inside of the list passed as last_row to go into the next dimension.
            store(array, indices, eq, dimension, iterations, new_row, interval_width) # Recursively call the store() function, this time with new_row in place of last_row.
            last_row.append(new_row) # Add new_row to last_row.
        else:
            # Unpack the indices using the * operator
            last_row.append(eq(*indices).integ(interval_width)) # If iterations == dimension, we want to evaluate our Galerkin Integral using its expression (eq()), our PofY objects (indices), and the integ() method with our interval_width.

        indices.pop() # Remove the last element from indices so that we can go through the for loop again with a different i value as the last element.



def solve_equation(equation_str, num_solutions, initial_guess, increments, equation_str2, resolution, num_variables, equation_str3):



    # This initializes the symbols for x and y
    x, y = symbols('x y')

    # equation(y, equation_str) just returns a lambdified equation for equation_str. The definition for this function is defined on line 65.
    eq = equation(y, equation_str)

    roots = []  # Initialize an empty list to store roots (the solutions for the characteristic equation)

    for i in range(num_solutions):
        # Use fsolve() to find the roots
        root = fsolve(eq, initial_guess, maxfev=1000000)[0]  # fsolve returns a 1-element list
        roots.append(root) # Add solution to the roots list.

        initial_guess += increments # increment the initial guess by the user-specified increments.

    print(f"{num_solutions} solutions: {roots}") # Print solutions.



    # What will happen is the computer will take the y values from the roots array and sub them in for y in the mode shape. The x values will be numbers ranging from 0 to 1 inclusive.
    # The number of x values is the resolution inputted by the user. These values are stored in the x_span array.
    x_span = np.linspace(0,1,int(resolution))

    # The following 5 lines are the creation of lambdified functions for the mode shape as well as its first, second, third, and forth derivatives.
    eq2 = lambdify((x, y), equation_str2, 'numpy')
    deq2 = lambdify((x, y), diff(equation_str2, x), 'numpy')
    ddeq2 = lambdify((x, y), diff(diff(equation_str2, x), x), 'numpy')
    dddeq2 = lambdify((x, y), diff(diff(diff(equation_str2, x), x), x), 'numpy')
    ddddeq2 = lambdify((x, y), diff(diff(diff(diff(equation_str2, x), x), x), x), 'numpy')

    # This is where we create lists to store the solutions to the lambdified equations above. These will be lists of lists, making them behave similar to 2D arrays.
    # The rows of these 2D arrays will correspond to each y value in the roots list and the columns will correspond to x values in x_span.
    matrix = []
    dmatrix = []
    ddmatrix = []
    dddmatrix = []
    ddddmatrix = []

    # The matrixify() function is used to calculate and store the solutions to the mode shape and it's derivatives in their respective arrays.
    # The function parameters are: The matrix we want to store the values in; the equation we are solving for; the roots; and x_span. You can find the function definition on line 68.
    print("Equation results: ")
    matrixify(matrix, eq2, roots, x_span)
    print("1st derivative results: ")
    matrixify(dmatrix, deq2, roots, x_span)
    print("2nd derivative results: ")
    matrixify(ddmatrix, ddeq2, roots, x_span)
    print("3rd derivative results: ")
    matrixify(dddmatrix, dddeq2, roots, x_span)
    print("4th derivative results: ")
    matrixify(ddddmatrix, ddddeq2, roots, x_span)

    # The rest of the solve_equation() function deals with creating a Galerkin Matrix. This matrix can have as many dimensions as the user wants.



    # This initializes the symbols for the Galerkin integral we will be making. The number of symbols is the same as the number of dimensions we want.
    P_symbols = symbols([f'P{j+1}' for j in range(num_variables)])

    eq = lambdify(P_symbols, equation_str3, "numpy") # Lambdify equation_str3 and store it in eq.

    # On lines 157-159, matrix2 is created. matrix2 is a list of objects of type PofY. PofY is the class defined at the beginning of the script.
    matrix2 = []
    for i, val in enumerate(matrix):
        matrix2.append(PofY(val, dmatrix[i], ddmatrix[i], dddmatrix[i], ddddmatrix[i]))

    interval_width = 1 / (len(x_span) - 1)
    K = [] # This is the list where the Galerkin matrix will be stored. This matrix can have more than 2 dimensions. The data will be stored in lists inside of lists inside of lists and so on, so that the K list can behave as if it has any number of dimensions.

    indices = [] # The indices list, in this case, will be different than the one used in the make_expr() function in the 2NSpace file. This one is going to store PofY objects.

    for j in matrix2:
        indices.append(j)
        iterations = 1 # Iterations is the number of dimensions K has so far. Keep in mind that when we pass it into the store() function on line 172, it is passed by value and is only changed INSIDE of the function.
        if iterations < num_variables: # num_variables is the same as the number of dimensions for the Galerkin matrix
            row = [] # If the number of dimensions currently in K is less than we want, we want to put a new list inside of K to give it another dimension.

            store(matrix2, indices, eq, num_variables, iterations, row, interval_width) # Call the store() function defined on line 76. Keep in mind that the row list is passed as the last_row parameter.
            K.append(row) # Append the row list to K.
        else:
            K.append(eq(*indices).integ(interval_width)) # If the Galerkin matrix only has one dimension, Then simply just find the integral for each element in matrix2 and append it to K.
        indices.pop() # Remove the last element from indices so that we can go through the for loop again with a different i value as the last element.
    print("\nK: ")
    print(K)
    return K, num_variables


if __name__ == "__main__":
    print("running\n")
    # This gets the user input for the characteristic equation. The equation is always 0 = (expression inputted by the user). Typically, this will be a trigonometric/hyperbolic-trigonometric equation.
    equation_str = input("Enter the expression for the characteristic equation in terms of y (use sin, cos, sinh, cosh): ")

    # Get user input for the number of solutions they want the computer to solve for the above equation.
    num_solutions = int(input("Enter the number of solutions you want: "))

    # Get user input for the initial guess since the computer will use the newton-raphson method to calculate solutions.
    initial_guess = float(input("Enter the initial guess: "))

    # This is the number of increments that the user wants to increment the initial guess by each time the computer calculates another solution.
    increments = float(input("Enter increments: "))

    # This is where the user inputs the expression for the mode shape in terms of x and y.
    equation_str2 = input("Enter right side of the mode shape in terms of x and y: ")

    # This is essentially the user telling the computer how many solutions for the mode shape they want.
    resolution = input("Enter desired resolution: ")

    # This is where the user inputs the number of dimensions they want in the Galerkin matrix.
    num_variables = int(input("Enter the number of desired dimensions for the Galerkin Matrix: "))


    # Here we ask the user to input the Galerkin integral.
    equation_str3 = input(f"Enter an expression in terms of {', '.join([f'P{j+1}' for j in range(num_variables)])}:" )

    # Call the solve_equation() function defined on line 92. This function returns the Galerkin matrix (K) along with its number of dimensions.
    K, dim = solve_equation(equation_str, num_solutions, initial_guess, increments, equation_str2, resolution, num_variables, equation_str3)
    print("Dimension: " + str(dim))

    # Printing the First Layer
    first_layer = 0
    for i in K:
        first_layer += 1
    print("First layer: " + str(first_layer))