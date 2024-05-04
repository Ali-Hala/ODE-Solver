import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy import symbols, lambdify, simplify

# Define the ODE system with user-specified number of variables
def model(t, y, ode_exprs):
    dydt = [ode_expr(y, t) for ode_expr in ode_exprs]
    return dydt

# Get user input for the number of variables
num_variables = int(input("Enter the number of variables: "))

# Get user input for the ODE expressions
ode_exprs = []
y_symbols = symbols([f'y{j+1}' for j in range(num_variables)])
t = symbols('t')

for i in range(num_variables):
    expr_str = input(f"Enter ODE expression for y{i+1} in terms of {', '.join([f'y{j+1}' for j in range(num_variables)])} and t: ")
    ode_expr = simplify(expr_str)
    ode_expr = lambdify((y_symbols, t), ode_expr, 'numpy')
    ode_exprs.append(ode_expr)

# Set initial conditions
y0 = [float(input(f"Enter initial condition for y{i+1}: ")) for i in range(num_variables)]

end = float(input("Enter the end:"))

# Get user input for the resolution (number of time points)
resolution = int(input("Enter the desired resolution (number of time points): "))

# Solve the ODEs using solve_ivp with the desired resolution
solution = solve_ivp(model, (0, end), y0, args=(ode_exprs,), method='RK45', t_eval=np.linspace(0, end, resolution), atol=1e-8, rtol=1e-8)

# Print the vectors corresponding to each variable and t
print("Time values (t):", solution.t)

for i in range(num_variables):
    print(f"Solution values (y{i+1}):", solution.y[i])

entry_count = 0
for i in solution.t:
    entry_count += 1

# Plot the results
for i in range(num_variables):
    plt.plot(solution.t, solution.y[i], label=f'y{i+1}')

print("Number of entries: " + str(entry_count))
plt.xlabel('Time')
plt.ylabel('y(t)')
plt.title('Numerical Integration of ODEs')
plt.legend()
plt.show()
