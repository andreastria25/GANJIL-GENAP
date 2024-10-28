import numpy as np
import matplotlib.pyplot as plt

# Diketahui
L = 0.5  # Henry
C = 10e-6  # Farad
fd = 1000  # Hz

# Fungsi untuk menghitung frekuensi f(R)
def f_R(R):
    return (1 / (2 * np.pi)) * np.sqrt((1 / (L * C)) - (R**2 / (4 * L**2)))

# Fungsi F(R) = f(R) - fd
def F_R(R):
    return f_R(R) - fd

# Check values of F_R at the endpoints
print("F_R(10):", F_R(10))
print("F_R(100):", F_R(100))

# Visualize F(R) to find a suitable interval
R_values = np.linspace(0, 200, 1000)
F_R_values = F_R(R_values)

plt.plot(R_values, F_R_values)
plt.axhline(0, color='red', lw=0.5)  # Horizontal line at y=0
plt.xlabel('R value')
plt.ylabel('F(R)')
plt.title('Plot of F(R)')
plt.grid(True)
plt.show()

# Implementasi metode biseksi
def bisection_method(func, a, b, tol=1e-6):
    if func(a) * func(b) > 0:
        raise ValueError("Fungsi tidak berubah tanda pada interval ini.")
    
    mid = (a + b) / 2.0
    R_values_bisection = []
    while (b - a) / 2.0 > tol:
        mid = (a + b) / 2.0
        R_values_bisection.append(mid)
        if func(mid) == 0:
            return mid, R_values_bisection
        elif func(a) * func(mid) < 0:
            b = mid
        else:
            a = mid
    return mid, R_values_bisection

# Mencari nilai R menggunakan metode biseksi
R_bisection, R_values_bisection = bisection_method(F_R, 0, 200)  # Adjust the interval

# Turunan dari F(R)
def dF_R(R):
    return -(R / (2 * np.pi * L**2 * np.sqrt((1 / (L * C)) - (R**2 / (4 * L**2)))))

def newton_raphson(func, dfunc, x0, tol=1e-6, max_iter=100):
    x = x0
    R_values_newton = []
    for i in range(max_iter):
        x_new = x - func(x) / dfunc(x)
        R_values_newton.append(x_new)
        if abs(x_new - x) < tol:
            return x_new, R_values_newton
        x = x_new
    raise ValueError("Metode Newton-Raphson gagal konvergen.")

# Mencari nilai R menggunakan metode Newton-Raphson
R_newton, R_values_newton = newton_raphson(F_R, dF_R, 50)

# Visualisasi konvergensi
def plot_comparison():
    plt.plot(R_values_bisection, label="Bisection Method")
    plt.plot(R_values_newton, label="Newton-Raphson Method")
    plt.xlabel('Iteration')
    plt.ylabel('R value')
    plt.title('Convergence Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()

plot_comparison()

# Fungsi eliminasi Gauss
def gauss_elimination(A, B):
    n = len(B)
    for i in range(n):
        for j in range(i+1, n):
            ratio = A[j][i] / A[i][i]
            for k in range(n):
                A[j][k] = A[j][k] - ratio * A[i][k]
            B[j] = B[j] - ratio * B[i]
    
    X = [0 for i in range(n)]
    X[n-1] = B[n-1] / A[n-1][n-1]
    
    for i in range(n-2, -1, -1):
        X[i] = B[i]
        for j in range(i+1, n):
            X[i] = X[i] - A[i][j] * X[j]
        X[i] = X[i] / A[i][i]
    
    return X

# Sistem persamaan
A = [[4, -1, -1], [-1, 3, -1], [-1, -1, 5]]
B = [5, 3, 4]

# Penyelesaian dengan eliminasi Gauss
solutions = gauss_elimination(A, B)
print(solutions)
