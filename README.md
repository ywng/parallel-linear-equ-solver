# parallel-linear-equation-solver
https://en.wikipedia.org/wiki/Jacobi_method

Let

{\displaystyle A\mathbf {x} =\mathbf {b} } A\mathbf {x} =\mathbf {b} 
be a square system of n linear equations, where:

{\displaystyle A={\begin{bmatrix}a_{11}&a_{12}&\cdots &a_{1n}\\a_{21}&a_{22}&\cdots &a_{2n}\\\vdots &\vdots &\ddots &\vdots \\a_{n1}&a_{n2}&\cdots &a_{nn}\end{bmatrix}},\qquad \mathbf {x} ={\begin{bmatrix}x_{1}\\x_{2}\\\vdots \\x_{n}\end{bmatrix}},\qquad \mathbf {b} ={\begin{bmatrix}b_{1}\\b_{2}\\\vdots \\b_{n}\end{bmatrix}}.} A={\begin{bmatrix}a_{11}&a_{12}&\cdots &a_{1n}\\a_{21}&a_{22}&\cdots &a_{2n}\\\vdots &\vdots &\ddots &\vdots \\a_{n1}&a_{n2}&\cdots &a_{nn}\end{bmatrix}},\qquad \mathbf {x} ={\begin{bmatrix}x_{1}\\x_{2}\\\vdots \\x_{n}\end{bmatrix}},\qquad \mathbf {b} ={\begin{bmatrix}b_{1}\\b_{2}\\\vdots \\b_{n}\end{bmatrix}}.


Input: initial guess {\displaystyle x^{(0)}} x^{{(0)}} to the solution, (diagonal dominant) matrix 
  
    
    {\displaystyle A}
  
A, right-hand side vector 
  
    
    {\displaystyle b}
  
b, convergence criterion
Output: solution when convergence is reached
Comments: pseudocode based on the element-based formula above
{\displaystyle k=0} k=0
while convergence not reached do
    for i := 1 step until n do
      {\displaystyle \sigma =0} \sigma =0  
      for j := 1 step until n do
        if j ≠ i then
          {\displaystyle \sigma =\sigma +a_{ij}x_{j}^{(k)}} \sigma =\sigma +a_{ij}x_{j}^{(k)}
        end
      end
      {\displaystyle x_{i}^{(k+1)}={{\frac {1}{a_{ii}}}\left({b_{i}-\sigma }\right)}} {\displaystyle x_{i}^{(k+1)}={{\frac {1}{a_{ii}}}\left({b_{i}-\sigma }\right)}}
    end
    {\displaystyle k=k+1} k=k+1
end
