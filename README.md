# Vanilla Option Pricing via Finite Difference Method

This project provides a tool to compute the price and Greeks of European and American vanilla options by solving the Black–Scholes PDE using a finite difference method (Crank–Nicolson scheme).  
The computational engine is written in C++, and the user interface is provided through an Excel (.xlsm) file.

## Features

- Price computation for European and American options (Call/Put)
- Greeks: Delta, Gamma, Theta, Rho, Vega
- Validation against the Black–Scholes closed-form solution
- Graphs:
  - Option Price vs Spot Price
  - Option Delta vs Spot Price
  - Exercise boundary for American options

## How to Use

1. **Download the files:**
   - The Excel interface file (`.xlsm`)
   - The dynamic library (`.dll`)

2. **Update the DLL path:**
   - Open the `.xlsm` file in Excel.
   - Open the VBA editor (shortcut: `Alt + F11`).
   - In *Module 1*, update the path to the downloaded DLL file.

3. **Input Parameters:**
   - Contract type (Call/Put)
   - Exercise type (European/American)
   - Maturity
   - Strike price
   - Computation date
   - Time mesh and spot mesh parameters
   - Current underlying price
   - Risk-free interest rate (as a piecewise linear function)
   - Volatility

4. **Run the Computation:**
   - Use the provided Excel interface to input the parameters and launch the pricing computation.

## Computation Method

- Solves the Black–Scholes PDE via a **finite difference method**.
- Crank–Nicolson time discretization.
- Special handling for American options (free boundary problem).

## Validation

- Results compared with the analytical Black–Scholes formula.
- Additional verification: Equivalence between European and American Call prices when \( r > 0 \) and Put prices when \( r < 0 \).

## References

- G. Allaire & S. M. Kaber, *Numerical Linear Algebra*, Springer, 2007
- S. Crépey, *Financial Modeling*, Springer, 2013
- G. Golub & C. Van Loan, *Matrix Computations*, Johns Hopkins University Press, 2013
- G. Pagès, *Numerical Probability*, Springer, 2018
- P. Wilmott, *Mathematics of Financial Derivatives: a Student Introduction*, Wiley, 1995
