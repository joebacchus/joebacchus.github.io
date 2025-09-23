---
title: Linear algebra, stochastic processes and optimisation
description: Interactive lecture notes for a third year mathematical methods course at Cambridge.
location: Cambridge
date: 2025-01-01
draft: true
---

# 3M1: Mathematical Methods - Linear Algebra, Stochastic Processes and Optimisation

Lent 2026

# Introduction

These notes cover four lectures and are just a taste of linear algebra. A number of topics that are important in scientific computing and engineering are touched upon, but there is much more left untold!

**Task**: Review carefully the linear algebra from Part IB.

## Books

There are many books on linear algebra. A gentle text is:

- Strang, G. (2006) *Linear Algebra and its Applications*.

An excellent textbook is

- Trefethen, L.N. and Bau, D (1997) *Numerical Linear Algebra*, SIAM.

Formal mathematical analysis of some of the topics covered can be found in:

- Süli, E. and Mayers, D. (2006) *An Introduction to Numerical Analysis*, Cambridge University Press.

Another book which also has a more formal emphasis, but which can be downloaded freely from the author's webpage is:

- Scott, L.R. (2014) *Numerical Analysis*, Princeton University Press. [http://people.cs.uchicago.edu/~ridg/newna/natwo.pdf](http://people.cs.uchicago.edu/~ridg/newna/natwo.pdf)

A classic matrix monograph is

- Golub, G.H. and Van Loan, C.F. (2012) *Matrix Computations*, The John Hopkins University Press.

# Definitions

Linear algebra involves a number of definitions and considerable jargon, and most problems and applications involve a synthesis of the basic concepts/definitions. Some key definitions and concepts that will be built upon in this course are presented in this section. Some concepts will be familiar, and others less so.

## Vector spaces

### What is a vector? [slide 2]

A vector $v$ is an element of a vector space. We will define a vector space shortly, but a space is not very intuitive without a concrete example of a vector.

You will be most familiar with vectors from the Euclidean space $\mathbb{R}^n$, which is all vectors with $n$ real components (an n-tuple). For example a vector $x \in \mathbb{R}^n$ has the form

$$ x = (x_1, x_2, x_3, \dots, x_n)^T $$

and a vector $x \in \mathbb{R}^4$ has the form

$$ x = (x_1, x_2, x_3, x_4)^T $$

For generality, we will consider complex problems. A vector $x \in \mathbb{C}^n$ has the form

$$ x = (x_1, x_2, x_3, \dots, x_n)^T $$

where $x_i \in \mathbb{C}$ (each component is a complex number).

### What is a vector space? [slide 3]

A vector space, often denoted by $V$, is a 'family' of vectors that obey some basic rules. For vectors $u, v$ and $w$ in a space $V$, and scalars $r$ and $s$ (complex), the rules are:

1. $u + (v + w) = (u + v) + w$
2. $u + v = v + u$
3. $u + 0 = u$
4. $u + (-u) = 0$
5. $r(s u) = (r s) u$
6. $r(u + v) = r u + r v$
7. $1 u = u$

Simply put, a vector in $V$ that is multiplied by a scalar is also in $V$. If the scalar is zero, the resulting vector is zero. Hence, a vector space must always contain a zero vector. The sum of any two vectors in vector space $V$ (possibly multiplied by a scalar) is also in the space.

### Vector space examples [slide 4]

**Example: real vectors in $n$ dimensions**
The vector space $\mathbb{R}^n$ is the space of all real vectors of length $n$.

**Example: polynomials of degree $p$**
All polynomials of degree $p$ in the variable $x$ on the interval $x \in (x_1, x_2)$ form a vector space. For example, if you add any two cubic polynomials, the results is also a cubic polynomial.

### What is a subspace? [slide 5]

A subspace of a vector space is a subset that obeys the rules of a vector space.

**Example: A subspace of $\mathbb{C}^3$**
All vectors of the form:

$$ x = \begin{pmatrix} x_1 \\ 0 \\ x_3 \end{pmatrix} $$

where $x_1$ and $x_3$ are any complex numbers, come from a subspace of $\mathbb{C}^3$. This subspace consists of all vectors that lie in the $x_1-x_3$ 'plane'.

Vectors of the form:

$$ x = \begin{pmatrix} x_1 \\ 1 \\ x_3 \end{pmatrix} $$

do not form a subspace since

$$ x + y = \begin{pmatrix} x_1 + y_1 \\ 2 \\ x_3 + y_3 \end{pmatrix} $$

Moreover, it does not contain the zero vector.

**Example: Algebra subspace of polynomials of degree $p$**
Polynomials of degree $p-1$ form a subspace of the space of all polynomials of degree $p$. A polynomial of degree $p$ contains all the lower order terms, i.e. a space of all quadratic polynomials contains all linear polynomials. The sum of two linear polynomials is also a linear polynomial.

## Matrices

### What is a matrix? [slide 6]

A matrix $A$ is a linear operator that performs a linear transformation from one vector to another:
$ A x = b $.
The above maps the vector $x \in \mathbb{C}^n$ to the vector $b \in \mathbb{C}^m$. We could also express this as $A: \mathbb{C}^n -> \mathbb{C}^m$. This says that the $A$ matrix maps a vector of length $n$ to a vector of length $m$.
A matrix is a rectangular array of numbers. For example, a $2 \times 3$ matrix has the form:

$$ A = \begin{pmatrix} A_{11} & A_{12} & A_{13} \\ A_{21} & A_{22} & A_{23} \end{pmatrix} $$

An $m \times n$ matrix comes from the space $\mathbb{C}^{m \times n}$, and we will often write $A \in \mathbb{C}^{m \times n}$.
A number of operations are defined for matrices, such as addition and multiplication, which you should be familiar with. One to should memorise is $C = A B$ using index notation:

$$ C_{ij} = \sum_k A_{ik} B_{kj} $$

Treating a column vector as an $n \times 1$ matrix, we can use the above to compute $A x = b$:

$$ b_i = \sum_j A_{ij} x_j $$

It is helpful to think about $b$ as a weighted sum of the columns of $A$, with the $j$th column multiplied by the $j$th entry in $x$.
Another perspective on matrix-matrix multiplication is via the sum of outer products,

$$ C = \sum_i a_i b_i^T $$

where $a_i$ is the $i$th column of $A$ and $b_i$ is the $i$th row of $B$. Note that $a_i b_i^T$ is a rank-1 matrix.

### What is the transpose and conjugate transpose? [slide 7]

The transpose of a matrix will be familiar,

$$ (A^T)_{ij} = A_{ji} $$

i.e. rows and columns are exchanged. In the context of complex problems, the conjugate transpose, also called hermitian transpose is more relevant:

$$ (A^H)_{ij} = \bar{A_{ji}} $$

where $\bar{x} \in \mathbb{C}$ is the complex conjugate of $x$. Expressed alternatively:

$$ A^H = \bar{A}^T = \bar{A^T} $$

For example,

$$ \begin{pmatrix} A_{11} & A_{12} & A_{13} \\ A_{21} & A_{22} & A_{23} \end{pmatrix}^H = \begin{pmatrix} \bar{A_{11}} & \bar{A_{21}} \\ \bar{A_{12}} & \bar{A_{22}} \\ \bar{A_{13}} & \bar{A_{23}} \end{pmatrix} $$

Clearly the transpose and the conjugate transpose coincide for real-valued matrices ($A \in \mathbb{R}^{m \times n}$).
A matrix $M \in \mathbb{C}^{n \times n}$ is Hermitian if

$$ M^H = M $$

if the matrix is real-valued, this is the same as the matrix being symmetric. Note that $(A B)^H = B^H A^H$. By extension, $A^H A$ must be Hermitian.
For real-valued matrices, the distinction between 'H' and 'T' vanishes.

### What is the dot/scalar/inner product? [slide 8]

The dot (or scalar) product is an operation between two equal-length vectors that yields a scalar, and is defined by:

$$ x^H y = \sum_{i=1}^n \bar{x}_i y_i $$

The dot product yields a complex number, and $x^H y = \bar{y^H x}$.
The notation $x \cdot y$ is sometimes used. An ambiguity is that this is sometimes defined to be

$$ x \cdot y = y^H x = \sum_{i=1}^n x_i \bar{y}_i $$

and other times

$$ x \cdot y = x^H y = \bar{y^H x} = \sum_{i=1}^n \bar{x}_i y_i $$

Since $x^H y \neq y^H x$ in the complex case and due to the ambiguity in the definition, we will not use the $\cdot$ notation.
The dot/scalar product is sometimes called the inner product as it is a inner product, but there are others operations that satisfy the requirements of an inner product.

*Complete Examples Paper questions 1 and 2.*

### What is an eigenpair? [slide 9]

For an $n \times n$ matrix $A$, $(\lambda, x)$ is an eigenpair of $A$ if

$$ A x = \lambda x $$

where $\lambda$ is an eigenvalue of $A$, and $x$ is the corresponding eigenvector of $A$. Recall that $\lambda$ can be equal to zero, but $x$ must be nonzero.
Recall that the eigenvalues of $A^{-1}$ are the reciprocal of the eigenvalues of $A$.

**Note**
Revise computation of eigenpairs for small matrices from Part I.

### Eigenvalues of a Hermitian matrix [slide 10]

Hermitian matrices have important properties that generalise those of symmetric matrices of real numbers:

- The eigenvalues of a Hermitian matrix are *real*.
- The eigenvectors of a Hermitian matrix are orthogonal, i.e. $u_i^H u_j = u_j^H u_i = 0$ when $i \neq j$.

(Proof in examples paper Q4)

### What is a unitary matrix? [slide 11]

Unitary matrices generalise the concept of orthogonal matrices of real numbers.
A matrix $Q \in \mathbb{C}^{n \times n}$ is a unitary matrix if $Q^H = Q^{-1}$, i.e. $Q^H Q = Q Q^H = I$.
A unitary matrix preserves the Euclidian norm of a vector.

### What is a Hermitian positive-definite matrix? [slide 12]

A Hermitian matrix $M \in \mathbb{C}^{n \times n}$ is positive definite if

$$ x^H M x > 0 \forall x \in \mathbb{C}^n \setminus \{0\} $$

The eigenvalues of a Hermitian positive definite matrix are strictly positive (this is a sufficient condition). Positive-definite matrices are particularly important for quadratic problems, which we will see later.

A Hermitian matrix $M \in \mathbb{C}^{n \times n}$ is semi-positive definite if

$$ x^H M x \geq 0 \forall x \in \mathbb{C}^n $$

which implies that all eigenvalues are positive.
The matrix $A^H A$ is positive semi-definite (see Examples Paper).

### What is the rank of a matrix? [slide 13]

The rank of a matrix $A \in \mathbb{C}^{m \times n}$ is the number of linearly independent rows or columns (the number is equal). It satisfies:

$$ \text{rank} A \leq \min(m, n) $$

A matrix is *full rank* if

$$ \text{rank} A = \min(m, n) $$

and *rank deficient* if

$$ \text{rank} A < \min(m, n) $$

### What is a sparse matrix? [slide 14]

A matrix in which most entries are zero is a sparse matrix. Typically, the number of non-zeroes on each row of a sparse matrix will be roughly the same and independent of the matrix size. Sparse matrices are very common in science and engineering, and often arise in connection with methods for solving differential equations, such as finite difference and finite element methods.

The number of entries in a dense $n \times n$ matrix is $n^2$. For a sparse matrix, the number of entries is $c n$, where $c$ is a small constant (say $c < 30$). For large $n$ the difference in the required storage is substantial.

**Exercise**
Compare the algorithmic complexity of matrix-vector multiplication for dense and sparse matrices.

*Complete Examples Paper questions 3, 4, 5 and 6.*

## Norms

### What is a norm? [slide 15]

A very important concept in linear algebra (and more generally vector spaces) is that of a norm. A norm of an object is a non-negative, real-valued number that is a measure of 'how big' something is and which allows ordering.
Norms that you will already be familiar with are the absolute value of a number, the modulus of a complex number and the Euclidean norm for vectors $x \in \mathbb{R}^n$:

$$ ||x|| = \sqrt{x_1^2 + x_2^2 + \dots + x_n^2} = (\sum_{i=1}^n x_i^2)^{1/2} $$

Denoting a norm of a vector $x$ by $||x||$ (not necessarily the Euclidean norm), a norm is a scalar that satisfies:

- $||x|| > 0$ when $x \neq 0$,
- $||x|| = 0$ when $x = 0$,
- $||k x|| = |k| ||x||$,
- $||x + y|| \leq ||x|| + ||y||$ (triangle inequality).

### Vector norms [slide 16]

A particular family of norms are known as $l_p$-norms:

$$ ||x||_p = (\sum_{i=1}^n |x_i|^p)^{1/p} $$

(recalling that $|x| = \sqrt{(\text{Re}(x))^2 + (\text{Im}(x))^2}$)). Commonly considered norms are the $l_1$ norm:

$$ ||x||_1 = |x_1| + |x_2| + \dots + |x_n| = \sum_{i=1}^n |x_i| $$

the $l_2$ norm (which we have already seen):

$$ ||x||_2 = (|x_1|^2 + |x_2|^2 + \dots + |x_n|^2)^{1/2} = (\sum_{i=1}^n |x_i|^2)^{1/2} $$

and the $l_\infty$ norm:

$$ ||x||_\infty = \max_i |x_i| $$

which is also known as the maximum norm.
The $l_1$ norm is sometimes called the 'taxicab norm'. Can you see why?
The Euclidean/$l_2$ norm is written as $||x||_2$, but since it is so frequently used the subscript '2' is sometimes dropped.
We can define norms that involve a matrix $A$, subject to some restrictions on the matrix:

$$ ||x||_A^2 = (x, A x) $$

where $(., .)$ is an inner product. Defining an inner product is not necessary at this stage, other than to say a common case is $(u, v) = v^H u$. We will see a common example shortly.

**Example: computing different norms for the same vector**
Compute the $l_1, l_2$ and $l_\infty$ norms for the vector

$ x = (2, -3, 7, 4) $

- $l_1$ norm: $||x||_1 = |2| + |-3| + |7| + |4| = 16$
- $l_2$ norm: $||x||_2 = (2^2 + 3^2 + 7^2 + 4^2)^{1/2} = \sqrt{78} \approx 8.83$
- $l_\infty$ norm: $||x||_\infty = \max_i |x_i| = 7$

### Example: different norms measure differently [slide 17]

Consider the two vectors

$ x = (3, 3, 3) $
$ y = (0, 0, 9) $

- Consider first the $l_1$ norm:
  $||x||_1 = |3| + |3| + |3| = 9$
  $||y||_1 = |9| = 9$
  In the $l_1$ norm $y$ and $x$ have the 'same magnitude'.

- Now the $l_2$ norm:
  $||x||_2 = (3^2 + 3^2 + 3^2)^{1/2} = \sqrt{27} \approx 5.20$
  $||y||_2 = (9^2)^{1/2} = \sqrt{81} = 9$
  In the $l_2$ norm $y$ is 'bigger' than $x$.

- Now the $l_\infty$ norm:
  $||x||_\infty = 3$
  $||y||_\infty = 9$
  In the $l_\infty$ norm $y$ is three times 'bigger than' $x$.

Note that $||x||_2^2 = x^H x$.

### Unit sphere and unit shell [slide 18]

It common to refer to to the *unit ball* or *unit sphere* for a given norm.

Unit ball (open)
$ \{ x \in V : ||x|| < 1 \} $

Unit ball (closed)
$ \{ x \in V : ||x|| \leq 1 \} $

Unit sphere
$ \{ x \in V : ||x|| = 1 \} $

The 'shape' depends on the chosen norm.

### Example: norm induced by a matrix [slide 19]

Linear systems $A x = b$ where the components of $x$ have units of length and the components of $b$ have units of force arise often in engineering and physics. It then follows that $x^H A x$ will have units of energy. We can use this to define a norm of $x$:

$$ ||x||_A^2 = x^H A x $$

For the above to obey the rules at the start of this section to qualify as a norm, $A$ must be positive definite. $A$ being positive definite implies that: (a) the energy is non-negative; and (b) the energy is zero only if $x = 0$. The above norm is often called the *energy norm*.

### Matrix norms [slide 20]

Norms can be defined for matrices, but they are less intuitive (at first) than norms of vectors, and can be more expensive to compute.

### Operator norms [slide 21]

A norm of a matrix $A$ is defined as:

$$ ||A|| = \max_{x \in \mathbb{C}^n \setminus \{0\}} (||A x||) / (||x||) $$

This norm measures the 'maximum amount' by which the matrix $A$ can re-scale a vector $x$ (in relative terms). From eq. (3), we can write:

$$ ||A x|| \leq ||A|| ||x|| \forall x $$

Matrix norms obey the rules at the start of section, and it follows from the definition that

$$ ||A B|| \leq ||A|| ||B|| $$

To quantify the 'size' of the original vector $x$ and the transformed vector $A x$, we need to choose a norm for the vectors.
Starting with $||A||_1$ (1-norm):

$$ ||A||_1 = \max_{x \in \mathbb{C}^n \setminus \{0\}} (||A x||_1) / (||x||_1) = \max_j \sum_{i=1}^n |a_{ij}| $$

which is the $l_1$-norm of the column of $A$ with the maximum $l_1$-norm.
For $||A||_\infty$ (infinity-norm):

$$ ||A||_\infty = \max_{x \in \mathbb{C}^n \setminus \{0\}} (||A x||_\infty) / (||x||_\infty) = \max_i \sum_{j=1}^n |a_{ij}| $$

which is the $l_1$-norm of the row of $A$ with the maximum $l_1$-norm. Proving the above two expressions is a question in the Examples Paper.
For $||A||_2$ (2-norm), we will first square both sides of eq. (3):

$$ ||A||_2^2 = \max_{x \in \mathbb{C}^n \setminus \{0\}} (||A x||_2^2) / (||x||_2^2) = \max (x^H A^H A x) / (x^H x) = \lambda_{\max}(A^H A) $$

where $\lambda_{\max}(A^H A)$ is the largest eigenvalue of $A^H A$ (recall that $A^H A$ is positive semi-definite, so all eigenvalues are positive).
The norm $||A||_2$ is therefore the square root of the largest eigenvalue of $A^H A$, or the largest singular values of $A$.
In the case that $A$ is Hermitian, the eigenvalues of $A^H A$ are the square of the eigenvalues of $A$. Hence, if $A$ is Hermitian $||A||_2 = |\lambda|_{\max}(A)$.

### Invariance of the 2-norm under rotation [slide 22]

The vector and matrix 2-norms are invariant under 'rotation', i.e. for $x \in \mathbb{C}^n$

$$ ||Q x||_2 = ||x||_2 $$

and for $A \in \mathbb{C}^{n \times m}$

$$ ||Q A||_2 = ||A||_2 $$

where $Q$ is a unitary matrix.

### Frobenius norm [slide 23]

Some matrix norms treat the entries of a matrix l(i k)e the entries of a vector. One such norm is the Frobenius norm. It is defined by:

$$ ||A||_F = \sqrt{\sum_i \sum_j |A_{ij}|^2} $$

The Frobenius norm is also invariant under rotation,

$$ ||Q A||_F = ||A||_F $$

where $Q$ is a unitary matrix.
Unless otherwise stated, when referring to matrix norms we mean operator norms.

### Which norm to choose? [slide 24]

Choosing a norm can depend on what we want to measure (what we want to 'weight' as important), and what fits naturally with a particular algorithm. We will see that some norms are easier to work with than others, so the choice is often pragmatic.
A technical point is that on the vector spaces that we consider, all norm are *equivalent*, which means that there exist constants $c_1 > 0$ and $c_2 > 0$ (but which typically depend on the dimension $n$) such that:

$$ c_1 ||x||_a \leq ||x||_b \leq c_2 ||x||_a \forall x \in V $$

# Stability and condition number

## Stability of operations

### Stability [slide 25]

Large linear systems of the form $A x = b$ are solved by computers. When solving with a computer, round-off error cannot be avoided. The important question is whether or not round-off errors will have a significant impact on the accuracy of the computed solution.
We can bound the error in terms of a quantity called the condition number.
The condition number of an invertible matrix $A$ is:

$$ \kappa(A) = ||A|| ||A^{-1}|| $$

To be concrete, we need to select a norm.

### Condition number in the 2-norm [slide 26]

For the 2-norm, we see that

$$ \kappa_2(A) = \sqrt{\lambda_{\max}(A^H A) / \lambda_{\min}(A^H A)} $$

since the eigenvalues of $A^{-1}$ are the reciprocal of the eigenvalues of $A$. if $A$ is Hermitian,

$$ \kappa_2(A) = |\lambda_{\max}(A)| / |\lambda_{\min}(A)| $$

## Error bounds

### Bounding the error when solving linear systems [slide 27]

Consider the problem

$$ A(x + \delta x) = b + \delta b $$

where $\delta b$ is the error in the RHS and $\delta x$ is the consequent error in the solution. Since $b = A x$ and $\delta x = A^{-1} \delta b$, we have

$$ ||\delta x|| / ||x|| \leq \kappa(A) ||\delta b|| / ||b|| $$

This shows how an error in $b$ can propagate through to the solution $x$. For large condition numbers we can expect small errors in $b$ to cause large errors in $x$.
An alternative, and sometimes more relevant scenario, is an error in the matrix $A$:

$$ (A + \delta A)(x + \delta x) = b $$

For example, the term $\delta A$ could represent the floating point errors introduced when performing an LU decomposition. In this case,

$$ ||\delta x|| / ||x|| \leq \kappa(A) ||\delta A|| / ||A|| $$

A matrix with a large condition number is said to be *ill-conditioned*.
It is important to note the distinction between the determinant and the condition number of a matrix. A small determinant does not necessarily mean that a matrix is ill-conditioned, and a moderate determinant does not mean that a matrix well conditioned.
Some examples are presented below, and larger examples can be found in the online Jupyter notebooks, and in particular for the notoriously ill-conditioned Hilbert matrix.

### Example: ill-conditioned system [slide 28]

Consider the problem

$$ \begin{pmatrix} 1 & 2 \\ 2 & 4.0001 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} b_1 \\ b_2 \end{pmatrix} $$

- When $b = (2, 4.0001)^T$, the solution is $x = (0, 1)^T$.
- When $b = (2, 4)^T$, the solution is $x = (2, 0)^T$.

We see here that a small change in the RHS leads to a very significant change in the solution. The eigenvalues are 5 and $2 \cdot 10^{-5}$, and the condition number for the 2-norm is 250,000.

### Example: moderate determinant, large condition number [slide 29]

The $n \times n$ matrix
$ A = \begin{pmatrix} 1 & -1 & \dots \\ 0 & 1 & -1 & \dots \\ \dots & & & 1 \end{pmatrix} $ (upper triangular with 1 on diagonal, -1 on super-diagonal)
has $\det(A) = 1$ and $k_\infty(A) = n 2^{n-1}$, hence the condition number becomes very large as $n$ increases, despite the determinant remaining constant at one.

### Example: small determinant, small condition number [slide 30]

The very simple $n \times n$ diagonal matrix

$$ D = \text{diag}(0.1, 0.1, \dots, 0.1) $$

has $\det(D) = 10^{-n}$ and $k_D(D) = 1$.

**Conditioning and determinant**
Conditioning should not be confused with the determinant.

*Complete Examples Paper questions 7, 8, 9, 10 and 11.*

# Interpolation and least-squares methods

## Interpolation

Interpolation is fitting a function to a data set that passes through the data points.

### Polynomial interpolation [slide 31]

if we have $n$ data points in a two-dimensional space $(x, y)$ we can usually fit a polynomial with $n$ coefficients. Say we have 6 data points $f = (f_0(x_0, y_0), \dots, f_5(x_5, y_5))^T$. We can (hopefully) interpolate the data points with a polynomial of the form:

$$ f(x, y) = c_0 + c_1 x + c_2 y + c_3 x y + c_4 x^2 + c_5 y^2 $$

We have six equations (one for each data point)

$$ f_i = c_0 + c_1 x_i + c_2 y_i + c_3 x_i y_i + c_4 x_i^2 + c_5 y_i^2 $$

and we can solve

$$ \begin{pmatrix} 1 & x_0 & y_0 & x_0 y_0 & x_0^2 & y_0^2 \\ \dots & \dots & \dots & \dots & \dots & \dots \\ 1 & x_5 & y_5 & x_5 y_5 & x_5^2 & y_5^2 \end{pmatrix} \begin{pmatrix} c_0 \\ \dots \\ c_5 \end{pmatrix} = \begin{pmatrix} f_0 \\ \dots \\ f_5 \end{pmatrix} $$

### Issues with polynomial interpolation [slide 32]

- As long as $A$ can be solved, i.e. the points do not all lie on a line, we can find the interpolating polynomial.
- The matrix $A$ is known as the Vandermonde matrix. It is a notoriously ill-conditioned matrix, and the condition number grows with increasing polynomial degree.
- Polynomial interpolation can be very sensitive to the choice of interpolation points.

### Interpolating a function [slide 33]

We consider the Runge function $f = 1/(1 + 25x^2)$ and interpolate it at five equally spaced points $x \in [-1, 1]$. We can fit a polynomial of degree four to the points:

[Figure: Points of the Runge function interpolated by a polynomial]
**Note the oscillations towards the end of the interval.**
Condition number $k_2 = 2.35 \times 10^1$.

### Interpolating a function: higher degree [slide 34]

Perhaps a interpolating polynomial of degree 12 will provide a better approximation:

[Figure: Points of the Runge function interpolated by a polynomial]
Condition number $k_2 \approx 1.23 \times 10^5$.

### Interpolating a function: higher degree issue [slide 35]

The amplitude of the oscillations near the ends of the domain are more severe at higher order. The is know as the Runge effect, and is particularly severe for high order interpolating polynomials with equispaced points.
A second, less obvious issue, is that the condition number of the Vandermonde matrix becomes very large as the polynomial degree is increased (see Jupyter notebooks), making numerical instability l(i k)ely.

### Orthogonal polynomials [slide 36]

A remarkably rich area that we will only briefly touch upon is orthogonal polynomials. So far we have considered the mononomial base $1, x, x^2, \dots, x^{n-1}$ and polynomials of the form

$$ f = c_{n-1} x^{n-1} + c_{n-2} x^{n-2} + \dots + c_1 x + c_0 $$

where we would pick (or solve for in the case of interpolation) the coefficients $\{c_i\}$.
There are alternatives to the monomial basis and which have a range of fascinating properties (far too many to cover).

### Legendre polynomials [slide 37]

We will consider Legendre polynomials on the internal $[-1, 1]$. There are various expressions computing the Legendre polynomials. One expression for $P_n$, is:

$$ (n+1)P_{n+1}(x) = (2n+1)x P_n(x) - n P_{n-1}(x) $$

where $P_0 = 1$ and $P_1 = x$. The special feature of Legendre polynomial is that:

$$ \int_{-1}^1 P_m(x) P_n(x) d\text{if} x = 0 \text{if} m \neq n $$

i.e. Legendre polynomials of different degree are orthogonal to each other.

### Legendre vs mononomial basis [slide 38]

The Legendre polynomials:

[Figure: Legendre polynomials of degree up to 9]
appear remarkably different from the mononimial basis:
[Figure: Monomial basis up to degree 9]

### Interpolation using Legendre polynomial [slide 39]

Legendre polynomials of degree up to and including $n$ span the same space as $1, x, x^2, \dots, x^n$, so we can express any polynomial of degree $n$ as

$$ f = \alpha_n P_n(x) + \alpha_{n-1} P_{n-1}(x) + \dots + \alpha_0 P_0(x) $$

To find the $\{\alpha_n\}$ coefficients we can construct a generalised Vandermonde matrix $A$ and solve $A \alpha = y_p$.

### 'Vandermode matrix' with Legendre polynomials [slide 40]

The matrix with Legendre polynomials is

$$ A = \begin{pmatrix} P_n(x_0) & P_{n-1}(x_0) & \dots & P_0(x_0) \\ \dots & \dots & \dots & \dots \\ P_n(x_n) & P_{n-1}(x_n) & \dots & P_0(x_n) \end{pmatrix} $$

In exact arithmetic both bases would compute the same polynomial.
Comparing the Vandermonde matrix on $[-1, 1]$ with equispaced points and $n=20$:
- Legendre basis: $k_2 = 2.72 \times 10^5$
- Monomial basis: $k_2 = 7.25 \times 10^8$
(see Jupyter notebook)

### How good is a polynomial interpolant? [slide 41]

For an function $f$, where $p_n$ is the best approximation of $f$ amongst all polynomials of degree $n$ and $p_i$ is the degree-n polynomial interpolant of $f$, then

$$ ||f - p_i|| \leq (\Lambda(T) + 1) ||f - p_n|| $$

where $\Lambda(T)$ is the Lebesgue constant. It depends on the choice of interpolation points $T$.
Beyond the scope of 3M1, but necessary for the above is the concept of a function norm. The simplest (and most commonly used) is the $L^2$-norm:

$$ ||f||_2 = (\int_\Omega f^2 d\text{if} x)^{1/2} $$

It is clear from earlier examples that exhibited the Runge effect that a polynomial that interpolated a function can be a very poor fit to the function at other points.

### Interpolation with non-equispaced points [slide 42]

It turns out that we can mitigate the Runge effect by using non-equispaced sampling points, and particularly good points are the roots of some orthogonal polynomials. Note from the plot of the Legendre polynomials that the roots tend to cluster close to the ends of the interval.

### Interpolation with non-equispaced points [slide 43]

Interpolation of the Runge function with a polynomial of degree 25 and using (i) equispaced points and (ii) the roots of the Legendre polynomials as sampling points:

[Figure: Interpolation with Legendre roots vs equispaced points]

For equispaced points, the Lebesgue constant grows exponentially in $n$. For Chebyshev points (roots of a family of orthogonal polynomials) the Lebesgue constant grows logarithmically in $n$.

## Chebyshev polynomials

### Chebyshev polynomials: definition [slide 44]

To be used later, we introduce the Chebyshev polynomials:

$$ T_n(\cos \theta) = \cos(n \theta), \quad \theta \in [0, \pi] $$

or

$$ T_n(x) = \cos(n \arccos x), \quad x \in [-1, 1] $$

The Chebyshev polynomial $T_n$ has $n+1$ extreme values, and the values have alternating sign:
$||T_n||_\infty = 1$, $T_n(x_k) = (-1)^k$ where $x_k = \cos(\pi k / n)$, $k=0, \dots, n$.
and $n$ distinct roots in $[-1, 1]$:
$T_n(t_k) = 0$ where $t_k = \cos((2k-1)\pi/(2n))$, $k=1, \dots, n$.

### Chebyshev polynomials: orthogonality [slide 45]

Consider

$$ \int_0^\pi \cos(n \theta) \cos(m \theta) d \theta = \begin{cases} 0 & \text{if } n \neq m \\ \pi & \text{if } n=m=0 \\ \pi/2 & \text{if } n=m \neq 0 \end{cases} $$

Introducing the substitution $x := \cos \theta -> d\text{if} x = -\sin \theta d \theta = -\sqrt{1-x^2} d \theta$, we have:

$$ \int_{-1}^1 (T_n(x) T_m(x))/\sqrt{1-x^2} d\text{if} x = \begin{cases} 0 & \text{if } n \neq m \\ \pi & \text{if } n=m=0 \\ \pi/2 & \text{if } n=m \neq 0 \end{cases} $$

### Chebyshev polynomials are polynomials [slide 46]

Straightforwardly,
$ T_0(x) = 1 $
$ T_1(x) = x $
Using the identity $\cos(a \pm b) = \cos a \cos b \mp \sin a \sin b$, we can derive
$ \cos(n \theta) = 2 \cos \theta \cos((n-1)\theta) - \cos((n-2)\theta) $
leading to, for $n \geq 2$,
$ T_n(x) = 2x T_{n-1}(x) - T_{n-2}(x) $

### Chebyshev basis [slide 47]

[Figure: Chebyshev polynomials]

### A Chebyshev 'minmax' property [slide 48]

For $y \notin (-1, 1)$, the solution to
$ \min \{ ||p(x)||_\infty : p \in P_n, -1 \leq x \leq 1, p(y) = 1 \} $
is $1/|T_n(y)|$.

**Message**: on the interval $[-1, 1]$, of all polynomials of degree $n$ that are equal to 1 at the fixed point $y$ outside of interval, the polynomial $T_n(x)/T_n(y)$ has the smallest deviation from zero.

## Data fitting

### Over-determined systems: residual [slide 49]

To find a solution to the problem $A x = b$, the vector $b$ must lie in the column space of $A$ (look back over your IB notes if you need to review this concept). if $A$ is an $m \times n$ matrix and $m > n$ (a skinny matrix), we have more equations than unknowns and in general there will be no solution.
For some vector $x$, we can define a 'residual' vector $r$:

$$ r = A x - b $$

### Over-determined systems: norm minimisation [slide 50]

To find an approximate solution to $A x = b$, we can search for a vector $x$ that minimises the residual in some norm:

$$ \min_x ||r(x)|| = \min_x ||A x - b|| $$

This says: find $\hat{x}$ that minimises the residual $r$ in the chosen norm.
The residual is a vector, so we need a way to compare different residual vectors to decide if one is smaller than another. For this we need to pick a norm in which to minimise the residual. Different norms will give different results, and it turns out that one norm will be much easier to work with than others.

### Minimising the error in the $l_2$ norm [slide 51]

if we use the $l_2$-norm for the problem in eq. (6), we seek:

$$ \min_x ||A x - b||_2 $$

Squaring $||A x - b||_2$ and expanding with indice:
$ r(x) = ||A x - b||_2^2 = (A x - b)^H (A x - b) = x^H A^H A x - x^H A^H b - b^H A x + b^H b $
which is the least-squares error. Hence, minimisation in the $l_2$-norm is known as the least-squares method.
Equation (8) is clearly quadratic in $x$. To find $\hat{x}$ that minimises $r$, we can differentiate $r$ and set the derivative equation to zero.
This gives the classic least-squares problem, also known as the normal equations:

$$ A^H A \hat{x} = A^H b $$

Note that a matrix of the form $A^H A$ is known as a normal matrix.
Re-arranging the least-squares problem:

$$ \hat{x} = (A^H A)^{-1} A^H b = A^+ b $$

### Pseudoinverse

The matrix $A^+ = (A^H A)^{-1} A^H$ is known as the pseudoinverse or the Moore-Penrose inverse. The inverse $(A^H A)^{-1}$ can be computed when $A$ is full rank, i.e. the columns of $A$ are linearly independent (see Examples Paper). In this section, $A$ has been a 'skinny matrix', i.e. $m > n$. if $A$ is 'fat matrix' ($m < n$) and full rank, then $A^+ = A^H (A A^H)^{-1}$.
There are cases where $A$ is not full rank, in which case the problem is under-determined and there are multiple solutions to the least-squares problem. We will look at strategies to handle this case once we have covered the singular value decomposition.

### Example: fitting a linear polynomial to a dataset in two dimensions [slide 52]

At the start of this section, we considered polynomial interpolation of six data points in a two-dimensional space. We now consider a least-squares fit with the linear polynomial

$$ f(x, y) = c_0 + c_1 x + c_2 y $$

We have six equations of the form

$$ f_i(x_i, y_i) = c_0 + c_1 x_i + c_2 y_i $$

and we can solve the normal equations.

### Least-squares fit of points on the Runge graph [slide 53]

We consider the problem we looked at earlier for polynomial interpolation. if we fit a 10th-order polynomial to 20 points:

[Figure: Least-squares fit of 20 points (10th order)]
The fit to the data looks reasonable. Fitting a 5rd-order polynomial:
[Figure: Least-squares fit of 20 points (5th order)]
we can see some difference between the fit (solid line) and the data points. Fitting a quadratic function:
[Figure: Least-squares fit of 20 points (quadratic)]
The fit is expectedly very poor.

### Solving least squares problems in practice [slide 54]

To solve the least squares problem $A^H A x = A^H b$, if $A$ is full rank we could apply LU factorisation to $A^H A$ (or more specially, Cholesky factorisation since $A^H A$ is Hermitian). A problem is that $A^H A$ is notoriously ill-conditioned, and the conditioning deteriorates as the matrix becomes larger.
Least-squares methods are usually solved using specialised methods that can manage the ill-conditioning of $A^H A$. QR factorisation (recall from Part IB) is better than LU, but can also suffer from accuracy/stability problems.

**Note**
In practice, do not write your own least-squares solver. Use a specialised function or library.

*Complete Examples Paper questions 12 and 13.*

### Minimising in other norms

Least-squares/2-norm minimisation is popular because the residual is quadratic, so when we take the derivative the minimisation problem becomes linear and is straightforward to solve. Also, being quadratic it will have a unique minimum (if $A$ is full rank).
We could pick another norm in which to minimise the error. The use of the 1-norm has become popular in recent years for some applications (e.g., compressed sensing). It is known as least absolute deviations. An issue with least-squares fitting is that it is sensitive to data outliers; the square of the 'error' amplifies the effect on the fitted function. Minimisation in the 1-norm is less sensitive to outliers.
The difficulty with norms other than the 2-norm is that the minimisation problem is non-linear, and therefore more challenging to solve.

# Iterative methods for linear systems

## Direct and iterative methods

We often want to solve a system of the form $A x = b$. Solution methods can be categorised as direct or iterative.

### Direct methods [slide 55]

A direct method computes a solution to $A x = b$ in a known/predictable number of operations. In the absence of round-off errors, the solution is exact. Solving a system via LU decomposition ($A x = L U x = b$) is an example of a direct method.
- For an dense $n \times n$ matrix, LU decomposition has complexity $O(n^3)$, which makes it very expensive for large matrices.
- For a sparse $n \times n$, the LU complexity is typically in the range $O(n^{3/2})$ and $O(n^2)$. This is clearly better than $O(n^3)$, but still expensive for large $n$.
- A more specialised point, direct methods cannot generally handle singular equations, i.e. matrices with a non-trivial nullspace.
- For reasonably conditioned systems, direct solvers can be robust, i.e. will reliably compute a solution.
- Direct solvers do not scale well on parallel computers since the substitution steps are inherently serial.

### Iterative methods [slide 56]

An iterative method seeks an approximate solution via a series of steps (iterations), and is usually terminated when the error/residual falls below a prescribed threshold (in the chosen norm!).
There are many iterative methods for solving problems in linear algebra. In general, these methods can sometimes be very fast, but can be slow or fail abjectly.
For useful iterative solvers:
- In some cases $O(n)$ schemes are possible (this is known as an optimal solver), and iterative solvers can be orders of magnitude faster than direct solvers.
- Use less memory than direct solvers.
- Some methods can solve singular problems.
- Iterations can be terminated early if high accuracy is not required.
- Iterative are often (but not always) suited to large, parallel computers.

## An iterative method for the largest eigenvalue and eigenvector

### Power iteration [slide 57]

A classic iterative method for finding the eigenvector associated with the largest absolute eigenvalue of a matrix is known as power iteration. Recall that a vector $x \in \mathbb{C}^n$ can be expressed in terms of the $n$ eigenvectors of an $n \times n$ matrix:

$$ x = \sum_{i=1}^n \alpha_i u_i $$

where $\alpha_i \in \mathbb{R}$. if we multiple $x$ repeatedly by $A$:

$$ A^k x = \sum_{i=1}^n \alpha_i \lambda_i^k u_i $$

we see that if the largest eigenvalue is distinct the resulting vector will be aligned with the eigenvector of the largest eigenvector (if $x$ has a component in the direction of the eigenvector).
The power iteration method is demonstrated in the Jupyter notebooks, where it is seen shown that convergence can be slow.
if we have an approximation of the eigenvector associated with $\lambda_{\max}$, a question is how can we find an approximation of $\lambda_{\max}$? In Part IA, you used the scaling from $A^k x$ to $A^{k+1} x$. However this can be very unreliable. There is a better way.

### Rayleigh quotient [slide 58]

Say we have an estimate of an eigenvector $x$ of the matrix $A$, in which case: $A x \approx \lambda x$.
To find an estimate of the corresponding eigenvalue $\lambda^*$, we could pose a minimisation problem in the 2-norm:

$$ \min_{\lambda^* \in \mathbb{C}} ||A x - \lambda^* x||_2 $$

This is minimised when

$$ \lambda^* = R(A, x) = (x^H A x) / (x^H x) $$

where $R$ is known as the Rayleigh quotient.
The Rayleigh quotient is often defined as being for Hermitian matrices only, in which case it must be real-valued. For Hermitian matrices it has a number of special properties, including second-order accuracy for estimating eigenvalues. That is, if the error in the eigenvector is $O(\epsilon)$, then the error in the eigenvalue estimated via the Rayleigh quotient is $O(\epsilon^2)$.

## Stationary methods for $A x = b$

### Family of stationary methods [slide 59]

We start with a family of simple methods for finding approximate solutions to $A x = b$. We decompose the matrix operator such that $A = N - P$. Rather than solving the exact problem, we solve

$$ N x = b + P x $$

iteratively. We compute an approximate solution $x_{k+1}$:

$$ N x_{k+1} = b + P x_k $$

where $x_k$ is the most recent known estimate of the solution, and the above equation is solved to compute the for the new estimate $x_{k+1}$. The process is then repeated to hopefully converge to the exact solution.
The key is to split $A$ such that equations of the form $N x = f$ are easy (inexpensive) to solve. Classic examples are:
- Richardson iteration: $N = I$.
- Jacobi method: $N = \text{diag}(A)$.
- Gauss-Seidel: $N = L(A)$ is the lower triangular part of $A$ (including the diagonal).

To understand whether or not we can expect these iterative methods to work, we need to examine what happens to the error as we iterate.

### Convergence [slide 60]

Defining the error at the $k$th iterate $e_k = x_{\text{exact}} - x_k$, from eq. (10) we have:

$$ N e_{k+1} = P e_k $$

and therefore

$$ e_{k+1} = N^{-1} P e_k $$

To converge, we need the error to approach zero as the number of iterations $k$ increases. We express the error vector $e_0$ as a linear combination of the eigenvectors of $N^{-1} P$:

$$ e_0 = c_1 u_1 + \dots + c_n u_n $$

and the error $e_k$ as

$$ e_k = (N^{-1} P)^k e_0 = c_1 \lambda_1^k u_1 + \dots + c_n \lambda_n^k u_n $$

The method will converge only if the absolute value of every eigenvalue is less than one.
The largest absolute eigenvalue of a matrix $A$ is often denoted by $\rho(A)$ and is known as the spectral radius.
The stationary methods based on splitting will converge if

$$ \rho(N^{-1} P) < 1 $$

The smaller $\rho < 1$, the faster the convergence. Examples of the Jacobi and Gauss-Seidel methods for a $50 \times 50$ matrix are presented in the Jupyter notebooks.

## Conjugate gradient method

The conjugate gradient (CG) method is remarkably simple but very powerful method. It is a Krylov subspace method, which is a more powerful family of the methods than the stationary methods in the previous section. The CG method is for Hermitian positive definite matrices only. There are other Krylov methods for other matrix types.
The CG method is sometimes presented as a direct method, but is applied in practice usually as an iterative method. The reasons will be clear when we consider some of its properties.

### Conjugate gradient as a direct method [slide 61]

We wish to solve $A x = b$, where $A$ is a $n \times n$ Hermitian positive-definite matrix. Consider that we have a set $P$ of $n$ non-zero vectors that are A-conjugate:

$$ P = \{p_0, p_1, \dots, p_{n-1}\} $$

A-conjugate implies that

$$ p_i^H A p_j = 0 \text{ if } i \neq j $$

Since $A$ is positive definite:

$$ p_i^H A p_i > 0 $$

Using $P$ as a basis for the solution:

$$ x = \sum_{i=0}^{n-1} \alpha_i p_i $$

and

$$ A x = \sum_{i=0}^{n-1} \alpha_i A p_i = b $$

Multiplying by $p_j^H$:

$$ p_j^H A x = \sum_{i=0}^{n-1} \alpha_i p_j^H A p_i = \alpha_j p_j^H A p_j = p_j^H b $$

Since the vectors in $P$ are A-conjugate:

$$ \alpha_j = (p_j^H b) / (p_j^H A p_j) $$

A question is how to generate the A-conjugate set $P$. A simple approach is to pick $n$ linearly-independent vectors and apply the Gram-Schmidt process to build $P$. However, this approach is not practical for large problems. Building and storing $P$ is too expensive.

### Conjugate gradient as an iterative method [slide 62]

There is a remarkable algorithm for the CG method in which elements of $P$ are computed with a very short recurrence relation. Moreover, in many cases we do not need all of $P$; sometimes we can get a very accurate solution with few elements of $P$.
Below is the CG algorithm. We will not prove it and you are not expected to memorise it, but it shows just how simple it is.

**Algorithm 1** Conjugate gradient method
1. $x_0 = 0, r_0 = b, p_0 = r_0$
2. for $k = 0, 1, \dots$ do
3. $\alpha_k = (r_k^H r_k) / (p_k^H A p_k)$
4. $x_{k+1} = x_k + \alpha_k p_k$
5. $r_{k+1} = r_k - \alpha_k A p_k$
6. $\beta_k = (r_{k+1}^H r_{k+1}) / (r_k^H r_k)$
7. $p_{k+1} = r_{k+1} + \beta_k p_k$
8. end for

- The only 'serious' operations you need to perform are matrix-vector products and dot products. For very large sparse problems, these two operations are amongst the easiest when using a parallel computer.
- Note also that only a few work vectors need to be stored. Contrast this with LU where the factorisations must be stored.
- We don't have time to derive the method, but we do want to know some of its key properties.

### Convergence [slide 63]

There is a lot of very rich analysis for the conjugate gradient method. We present here some key results.
Solving $A x = b$ where $A$ is an $n \times n$ matrix:
- For the error $e_k = x_{\text{exact}} - x_k$, the CG method is monotone in the A-norm: $||e_{k+1}||_A \leq ||e_k||_A$.
- From the preceding point, the CG method will solve the problem exactly (in the absence of round-off error) in at most $n$ iterations. This is why it is sometimes considered a direct method.
- The number of iterations required to solve exactly equal to the number of distinct eigenvalues of $A$.
- The rate of convergence is affected by the condition number $\kappa_2(A)$:

$$ ||e_k||_A / ||e_0||_A \leq 2 ((\sqrt{\kappa_2(A)} - 1)/(\sqrt{\kappa_2(A)} + 1))^k $$

**Example**
The monotonic convergence is demonstrated by the example in the Jupyter notebooks.

### Preconditioning [slide 64]

if the condition number of a matrix is large, the CG method may be too slow to converge for practical use. In this case, preconditioning can be attempted, where the transformed system:

$$ P^{-1} A x = P^{-1} b $$

This is known as 'left preconditioning'.

### Preconditioning idea [slide 65]

- if $P^{-1} \approx A^{-1}$, the condition number of the $P^{-1} A$ will be better than $A$, and the CG method will therefore converge faster.
- In practice it can be a balancing act - $P^{-1}$ must be cheap to apply for efficiency, but it must be close enough to $A^{-1}$ to be effective.
- When it all works, preconditioned CG solvers can be orders of magnitude faster than LU factorisations, and they work much better on large, parallel computers.

*Complete Examples Paper question 14.*

# Singular value decomposition (SVD)

Singular value decomposition (SVD) is one of the most beautiful concepts in mathematics. The definition is wonderfully simple, but the insights and applications are rich.

## Definition and properties

### Matrix diagonalisation – review [slide 66]

The diagonalisation of real, symmetric matrices was covered in Part IA. Extending this to Hermitian matrices, for a Hermitian matrix $M \in \mathbb{C}^{n \times n}$:

$$ Q^H M Q = \Lambda $$

where the columns of $Q$ are the normalised (in $l_2$) eigenvectors of $M$ and $\Lambda$ is a diagonal matrix of the eigenvalues of $M$ (which are real).
Since the eigenvectors of a Hermitian matrix are orthogonal, $Q$ is a unitary matrix, i.e. $Q^{-1} = Q^H$, and

$$ M = Q \Lambda Q^H $$

Matrix diagonalisation can be very helpful because we see that in a particular basis ('co-ordinate system') the action of a matrix on a vector is easy to interpret. Moreover, some operations, such as matrix-matrix multiplication, become trivial.
Whilst matrix diagonalisation is a very powerful concept it has some major limitations:
- It is only valid for square matrices;
- Not all matrices can be diagonalised, e.g. defective matrices; and
- The matrix $Q$ of normalised eigenevctors is guaranteed unitary only for Hermitian matrices.

**Multplication by a diaginal matrix**
if $D$ is a diagonal matrix wih diagonal entries $d_1, \dots, d_n$, then

$$ A D = (d_1 a_1, \dots, d_n a_n) $$

where $a_i$ is the $i$th column of $A$.

**Expanding the diagonalisation**
From the outer product representation of a matrix-matrix product and eq. (12), expanding eq. (11) we get

$$ M = \sum_i \lambda_i u_i u_i^H $$

where $(\lambda_i, u_i)$ is the $i$th eigenpair of $M$. The term $u_i u_i^H$ is a rank-1 matrix.

**Hermitian matrices**
All Hermitian matrices can be diagonalised

**Sign ambiguity**
Recall that eigenvectors of a matrix are not unique; it is direction that is important. Even when normalised, there persists a sign ambiguity:

$$ M = Q \Lambda Q^H = (-Q) \Lambda (-Q)^H $$

### Definition [slide 67]

The singular value decomposition of an $m \times n$ matrix $A$ is:

$$ A = U \Sigma V^H $$

where
- $U \in \mathbb{C}^{m \times m}$ is a unitary matrix;
- $\Sigma \in \mathbb{R}^{m \times n}$ is a diagonal matrix, with diagonal entries $\sigma_i$ (the 'singular values') sorted such that $\sigma_1 \geq \sigma_2 \geq \dots \geq \sigma_\rho \geq 0$, where $\rho = \min(m, n)$; and
- $V \in \mathbb{C}^{n \times n}$ is a unitary matrix.

if $A$ was Hermitian, then we would have $U = V = Q$, and $\Sigma = \Lambda$.

### Question: what should $U, \Sigma$ and $V$ be? [slide 68]

1. To answer part of this question, premultiply both sides of eq. (14) by $A^H$:

$$ A^H A = (U \Sigma V^H)^H U \Sigma V^H = V \Sigma^H U^H U \Sigma V^H = V (\Sigma^H \Sigma) V^H $$

Noting that $A^H A$ is Hermitian and that $\Sigma^H \Sigma$ is diagonal with entries $\sigma_1^2, \sigma_2^2, \dots, \sigma_n^2$, comparing this to the diagonalisation of a Hermition matrix in eq. (11), we see that this is just diagonalisation of the square matrix $A^H A$. Therefore:
- The columns of $V$ are the (normalised) eigenvectors of $A^H A$; and
- The diagonal entries of $\Sigma^H \Sigma$ are the eigenvalues of $A^H A$.

Now post-multiplying both sides of eq. (14) by $A^H$:

$$ A A^H = U \Sigma V^H (U \Sigma V^H)^H = U \Sigma V^H V \Sigma^H U^H = U (\Sigma \Sigma^H) U^H $$

We now see that:
- The columns of $U$ are the (normalised) eigenvectors of $A A^H$; and
- The diagonal entries of $\Sigma \Sigma^H$ are the eigenvalues of $A A^H$, which are the same as the eigenvalues of $A^H A$.

What is missing is the signs of the eigenvectors. We consider this next.
2. Post-multiplying both sides of $A = U \Sigma V^H$ by $V$,

$$ A V = U \Sigma $$

which can also be expressed as

$$ A v_i = \sigma_i u_i $$

where $v_i$ is the $i$th column of $V$ and where $u_i$ is the $i$th column of $U$.

### SVD definition: summary [slide 69]

- The columns of $U$ are the (normalised) eigenvectors of $A A^H$;
- The diagonal entries of $\Sigma$ are the square roots of the eigenvalues of $A^H A$ (or, equivalently $A A^H$); and
- The columns of $V$ are the (normalised) eigenvectors of $A^H A$;
- Use $A v_i = \sigma_i u_i$ to deduce the sign for $u_i$ given $v_i$ (or vice versa); and
- Every matrix has a SVD, which is unique up to the signs of the eigenvectors in $U$ and $V$.

The SVD of a Hermitian matrix is the usual eigen diagonalisation.

### The reduced SVD [slide 70]

Visualising the shape of the SVD on an $m \times n$ matrix with $m > n$:

[Figure: Visualising the shape of the SVD]

the terms below row $n$ in $\Sigma$ are always zero. This means that the last $m-n$ columns of $U$ make no contribution. In practice, the 'reduced' SVD, in which redundant entries are removed, is typically used:

[Figure: Visualising the shape of the reduced SVD]

### Computing the SVD: example [slide 71]

How to compute an SVD robustly could be an entire lecture course on its own. We will now compute an example 'by hand' by computing eigenpairs of $A A^H$ and $A^H A$. As is often the case in linear algebra, the method by which we compute problems by hand is not the way real problems should be computed. As touched upon with least-squares methods, the normal matrix $A A^H$ is notoriously ill-conditioned. Use built-in library functions to compute an SVD.

[NOTE: There appear to be inconsistencies in the OCR for the matrices on this page.]
[The following is a direct transcription of the provided text.]

Compute the SVD of the matrix $A = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 1 \end{pmatrix}$.
1. Compute eigenvalues/vectors of $A A^T$:
$A A^T = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$.
For this matrix, $\lambda_1 = 3$ and $u_1 = (1/\sqrt{2}) \begin{pmatrix} -1 \\ 1 \end{pmatrix}$, and $\lambda_2 = 1$ and $u_2 = (1/\sqrt{2}) \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.
Therefore $U = 1/\sqrt{2} \begin{pmatrix} -1 & 1 \\ 1 & 1 \end{pmatrix}$ and $\Sigma = \begin{pmatrix} \sqrt{3} & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$.

2. Compute eigenvalues/vectors of $A^T A$:
$A^T A = \begin{pmatrix} 1 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 1 \end{pmatrix}$.
For this matrix, $\lambda_1=3, u_1=(1/\sqrt{6})\begin{pmatrix} 1 \\ -2 \\ 1 \end{pmatrix}$, $\lambda_2=1, u_2=(1/\sqrt{2})\begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix}$, and $\lambda_3=0, u_3=(1/\sqrt{3})\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$.

3. The last column of $\Sigma$ is zero and therefore make no contribution. The reduced SVD is therefore:
$A = 1/\sqrt{2} \begin{pmatrix} -1 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} \sqrt{3} & 0 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1/\sqrt{6} & -2/\sqrt{6} & 1/\sqrt{6} \\ -1/\sqrt{2} & 0 & 1/\sqrt{2} \end{pmatrix}^T$

### Low-rank approximations [slide 72]

if we expand the SVD in eq. (14), we get

$$ A = \sum_{i=1}^r \sigma_i u_i v_i^H $$

where $r$ is the number of non-zero singular values, $u_i$ is the $i$th column of $U$ and $v_i$ is the $i$th column of $V$. The above is the expression of a matrix as a sum of rank-1 matrices.
A 'low rank' approximation of $A$ is

$$ A_k = \sum_{i=1}^k \sigma_i u_i v_i^H $$

where $k < r$. The rank of $A_k$ is $k$, which is less than the rank of $A$ ($A_k$ has fewer linearly independent rows/columns than $A$).
How good is $A_k$? It is possible to show that

$$ ||A - A_k||_F = \sqrt{\sigma_{k+1}^2 + \dots + \sigma_r^2} \leq ||A - B||_F $$

for all matrices $B$ of rank $k$ or less. This says that $A_k$ is a better approximation of $A$ than any other matrix of rank $k$, measured in the Frobenius norm.
Similarly in the 2-norm:

$$ ||A - A_k||_2 = \sigma_{k+1} \leq ||A - B||_2 $$

for all matrices $B$ of rank $k$ or less. The key result is that SVD can be used to construct optimal low-rank approximations of a matrix.

## Applications of the SVD

There are many applications of the SVD. Here we present just a few. All the presented examples can be found in the Jupyter notebooks.

### Finding patterns in data [slide 74]

We can use SVD to interpret large data sets. As an example, consider a $100 \times 200$ matrix with entries that are equal to 0 or 1. The entries have been arranged in a pattern which is shown graphically below:

[Figure: Matrix with a pattern]

Performing a SVD of the matrix, we find $\{\sigma_i\} = \{1.35 \times 10^2, 1.96 \times 10^1, 1.32 \times 10^1, 0, \dots, 0\}$.
This means that the matrix is rank 3, and can be represented as the sum of just three rank-1 matrices.
We now take the same problem, and add noise (between 0 and 0.1) to the 'white' background image:

[Figure: Matrix with noise]

if you look carefully you can see squares for each matrix entry (depending on print quality, you might want to look at this example on a screen).
if we perform an SVD of the matrix with noise, the largest singular values is roughly $1.29 \times 10^2$, and the smallest is roughly $1.2 \times 10^{-1}$. Reconstructing the image using only singular values that are larger than 1 (there are three of these):

[Figure: Matrix after SVD cleaning]

Visually, much of the noise has been removed.

### Image compression [slide 75]

Images can be represented as matrices, with one matrix entry for each pixel. The matrix entry is a number that represents an intensity. For RGB colour images, an image is represented by three such matrices - one for red, one for green and one for blue. In 24-bit colour, each colour at each pixel is represented with 8 bits, which means an integer between 0 and 255 for the colour intensity.
An 8 megapixel digital image has $2448 \times 3264$ pixels. Storing this as an RGB matrix would require roughly 24MB. To compress the image, we could try a low-rank approximation. We can construct an approximation by computing the SVD of the image and discarding the singular values that fall below a threshold.
Below is a photograph taken in the first lecture:

[Figure 1: Original image.]

We can represent the grey scale the image as a matrix with dimensions equal to the number of pixels and with float values between 0 and 1 corresponding the intensity (0=black, 1=white), or integers between 0 and 255 (0=black, 255=white).
**Singular values** Performing a SVD of the image, we plot the singular values:

[Figure 2: Singular values of the image.]

The ratio between the largest and smallest singular values is several orders of magnitude.
We now create a low-rank approximation of the original image by retaining only the larger singular values. if we retain the largest 10% of the singular values:

[Figure 3: Compressed (left) with largest 10% of the singular values and original (right) images.]

if we retain only the largest 2% of the singular values:

[Figure 4: Compressed (left) with largest 2% of the singular values and original (right) images.]

if we retain only the largest 0.5% of the singular values:

[Figure 5: Compressed (left) with largest 0.5% of the singular values and original (right) images.]

### Effective rank [slide 76]

When taking measurements, noise is often unavoidable and this can make it hard to detect (near) linear dependencies. Consider the matrix

$$ A = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 2 & 2 \\ 1 & 0 & 1 \end{pmatrix} $$

By inspection, this matrix is clearly rank 2. if we add noise to the matrix entries in the range $(0, 10^{-6})$:
it becomes full rank (rank 3). In practice we would work with large data sets, so would have little chance of detecting any near linear dependencies by inspection.
Performing an SVD on the above matrix the singular values are: 4.048, 0.7811 and $3.90 \times 10^{-7}$. The effective rank is the number of singular values that are greater than the noise level, i.e. in this case the effective rank is 2.

### Least-squares solutions [slide 77]

We'll now do some linear algebra gymnastics to use the SVD to better understand least squares problems, including the rank deficient case.
Recall that if the $m \times n$ matrix $A$, with $m > n$, is full rank, then the least squares solution to $A x = b$ is:

$$ \hat{x} = (A^H A)^{-1} A^H b = A^+ b $$

where $A^+ = (A^H A)^{-1} A^H$ is the psuedoinverse.

### Full rank case [slide 78]

if $A$ is an $m \times n$ matrix with $m > n$, we can partition the SVD as:

$$ A = (U_1, U_2) \begin{pmatrix} \Sigma_1 \\ 0 \end{pmatrix} V^H $$

Since matrix is full rank, $\Sigma_1$ will be $n \times n$, and $U_1$ is $m \times n$. Since the columns of $U$ are orthonormal, $U_1^H U_1 = I_{n \times n}$ and $U_2^H U_1 = 0_{(m-n) \times n}$. Recall that $||A x||_2 = ||U A x||_2$ when $U$ is unitary. We will use this property below.
A least squares solution minimises:

$ r = ||A x - b||_2^2 = ||U^H(A x - b)||_2^2 = ||\begin{pmatrix} \Sigma_1 V^H x - U_1^H b \\ -U_2^H b \end{pmatrix}||_2^2 = ||\Sigma_1 V^H x - U_1^H b||_2^2 + ||U_2^H b||_2^2 $

The residual is obviously minimised when $\Sigma_1 V^H x = U_1^H b$, i.e.

$$ \hat{x} = V \Sigma_1^{-1} U_1^H b $$

Both $\Sigma_1$ and $V$ are full rank, therefore the least-squares solution is unique.
if we wish, we can add the 'zero padding' back to $U$:

$$ \hat{x} = V \Sigma^+ U^H b $$

where $\Sigma^+$ is the pseudo inverse of $\Sigma$: $\Sigma^+ = \begin{pmatrix} 1/\sigma_1 & \dots \\ \dots & 1/\sigma_n & \dots \end{pmatrix}$.
The SVD can be used to compute least-squares solutions. Compared to Cholesky factorisation of the normal matrix $A^H A$ or QR factorisation, it is the most stable.

### Relationship to the normal equations [slide 79]

if $A$ is full rank and we have a SVD of $A$:

$$ A^+ = (A^H A)^{-1} A^H = (V \Sigma_1^2 V^H)^{-1} V \Sigma_1 U_1^H = V \Sigma_1^{-1} U_1^H $$

This gives the least squares solution $\hat{x} = V \Sigma_1^{-1} U_1^H b = V \Sigma^+ U^H b$.

### Stability [slide 80]

if $\hat{x} = A^+ b$ ($x \in \mathbb{C}^n$) then:

$$ ||\hat{x}||_2^2 = ||\Sigma_1^{-1} U_1^H b||_2^2 \geq |u_n^H b|^2 / \sigma_{\min}^2 $$

This says that if the smallest singular value $\sigma_{\min}$ is small, then the least squares solution will be large and very sensitive to changes in $b$.

### Rank deficient case [slide 81]

The preceding shows that if $\sigma_{\min} = 0$, we cannot compute a unique least squares solution. In the case that $\sigma_{\min}$ is small, while we can compute a solution, $||x||_2$ can become very large with small changes in $b$, which is not very satisfactory either.
For a problem with zero singular values, consider a partitioning of the SVD:

$$ A = (U_1, U_2) \begin{pmatrix} \Sigma_1 & 0 \\ 0 & 0 \end{pmatrix} (V_1, V_2)^H $$

where $\Sigma_1$ has size $r \times r$. It follows that

$$ A = U_1 \Sigma_1 V_1^H $$

Examining the least-squares residual (recall that we can apply $(U_1, U_2)^H$ to $A x - b$ and it will not change the $l_2$ norm):

$$ ||A x - b||_2^2 = ||\Sigma_1 V_1^H x - U_1^H b||_2^2 + ||U_2^H b||_2^2 $$

It is clear that $||A x - b||_2$ is minimised when $\Sigma_1 V_1^H x = U_1^H b$.
For the full rank case, we worked with $V^H$, which has only the trivial nullspace. However, we now work with $V_1^H$, which is only part of $V^H$, and the matrix $V_1^H$ consequently has a non-trivial nullspace. Specifically, we have $V_1^H(V_2 z) = 0$ for all $z \in \mathbb{C}^{n-r}$ since the columns of $V_1$ and $V_2$ are orthogonal. The result is that there are an infinite number of solutions $\hat{x}$ that satisfy $\Sigma_1 V_1^H x = U_1^H b$.
Since vectors $V_2 z$ lie in the nullspace of $V_1^H$, we have that

$$ \hat{x} = V_1 \Sigma_1^{-1} U_1^H b + V_2 z $$

for all $z \in \mathbb{C}^{n-r}$ is a solution to the least squares problem. Taking the norm of the above expression and exploiting that $V_1$ and $V_2$ are orthogonal with respect to each other:

$$ ||\hat{x}||_2^2 = ||V_1 \Sigma_1^{-1} U_1^H b||_2^2 + ||V_2 z||_2^2 $$

Therefore

$$ \hat{x} = V_1 \Sigma_1^{-1} U_1^H b $$

is the minimiser to the least-squares problem with minimal $l_2$ norm. In other words, of the all the solutions to the least squares problem, $\hat{x}$ is the smallest in the $l_2$ norm.

### Fitting points and singular systems [slide 82]

Say we are given four data points that depend on $x$ and $y$, and we are asked to fit a polynomial of the form

$$ f(x, y) = c_{00} + c_{10} x + c_{01} y + c_{11} x y $$

to the data points. Normally, we would expect to be able to fit the above polynomial to four data points by interpolation, i.e. solving $A c = f$ where $A$ a square Vandermonde matrix. However, if the points happen to lie on a line, then $A$ will be singular. if the points happen to lie almost on a line, then $A$ will be almost singular.
A possibility is to exclude zero or small singular values from the process, thereby finding a least-squares fit with minimal $||c||_2$. We test this for the data set $f_1(1,0)=3, f_2(2,0)=5, f_3(3,0)=7, f_4(4,0)=9$. The data lies on the line $y=0$, and is in fact is linear in $x$.
To find the polynomial coefficients we want to solve $A c = f$. For the points in the data set we have:

$$ A = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 1 & 2 & 0 & 0 \\ 1 & 3 & 0 & 0 \\ 1 & 4 & 0 & 0 \end{pmatrix} $$

which is clearly singular. We could try a least-squares fit by solving $A^T A c = A^T f$, but $A^T A$ is singular.
if we perform an SVD of $A$ and compute $c = V \Sigma^+ U^T f$, we get

$$ f(x, y) = 1 + 2x $$

which in fact interpolates the points.
See the Jupyter notebook for the steps.

### Example: fitting with nearly singular systems [slide 83]

if we take the previous example and add a small amount of noise $\epsilon \in (-5 \times 10^{-4}, +5 \times 10^{-4})$ to $x_i, y_i$ and $f_i$, the Vandermonde matrix is no longer singular and can be used to compute an interpolating polynomial. For some noise, we find that

$$ f(x, y) = 1.00365037 + 1.99853161x - 5.16321091y + 2.40001974 x y $$

The coefficients of the $y$ and $x y$ terms are suddenly significant, and this is due only to the noise.
if we compute and SVD and ignore singular values less than $10^{-3}$, we get:

$$ f(x, y) = 0.999257206 + 2.00013498x + 6.49882602 \times 10^{-4} y + 1.13588861 \times 10^{-3} x y $$

which is a more reasonable fit and very close to the noise-free case. See the Jupyter notebook for the steps.

*Complete Examples Paper questions 15, 16, 17 and 18.*

# Principal component analysis (PCA)

### What is principal component analysis (PCA)? [slide 84]

PCA is a dimension reduction technique.
- From a dataset, it attempts to find new, uncorrelated variables
- The 'new' variables are the directions along which the variance in the data is greatest
- Effectively, data is transformed/rotated (linearly) to a new coordinate systems, with the new coordinates ordered starting with the direction of greatest variation.
- Geometric interpretation: fits a (hyper-)ellipse to a dataset, the principal axes are the directions with the greatest variance (in descending order, starting from the longest).
- It finds linear subspaces – works well when 'good' new variables can be defined as linear combinations of the original variables.

### Example problem: arms and legs [slide 85]

Consider arm and leg length measurements for a group of $n$ people. We have $n$ samples and $p=2$ variates (arm length and leg length). Putting the data in a $n \times p$ matrix $X$.
Intuitively we would expect a strong correlation between arm and leg length. We can use PCA to uncover this.

### Arms and legs: coordinates [slide 86]

- We want to find directions in our 2D space (arm length, leg length) along which the variance is greatest. if there is one dominant direction we could reduce the number of parameters from 2 to 1.
- The component of the point $(X_{i1}, X_{i2})$ in the direction $a$, where $||a||_2 = 1$ is $d_i = (X_{i1}, X_{i2}) (a_1, a_2)^T$ (in the $l_2$-norm). if $d_i$ varies significantly as we move in the $(a_1, a_2)$ direction, then a coordinate in that direction is a good candidate for an independent variable in the problem.

### Arms and legs: perfect data [slide 87]

Imagine a dataset where $l_{\text{arm}}$ is the arm length and $l_{\text{leg}} = \alpha l_{\text{arm}}$, where $\alpha > 0$ is a constant. For each person we have the measurements $(l_{\text{arm}}, \alpha l_{\text{arm}})$. All points are spread along the line in the $(1, \alpha)$ direction.
In the orthogonal direction $(-\alpha, 1)$, there is no change in the data.
In practice data will have noise, correlations will no be perfect, and higher-dimensional problems cannot be easily visualised.

### Arms and legs: variance [slide 88]

Generalising $d_i = (X_{i1}, X_{i2}) (a_1, a_2)^T$, the component of each sample is the $a$ direction is:

$$ d = X a $$

The aim will be to find the direction $a$ that maximises the variance of $d$.
PCA finds the orthogonal directions $a_i$, from greatest to least variance.

### Covariance matrix [slide 89]

Let $X \in \mathbb{R}^p$ be a random vector. The (co)variance matrix is:

$$ \text{var}(X) := E((X - \mu)(X - \mu)^T) $$

The covariance matrix has shape $p \times p$. It is symmetric and positive-semidefinite.

### Sample data [slide 90]

Imagine that we have samples for $n$ entities/people, and for each entity/person we have $p$ variates (measurement variables), e.g. arm length and leg length for each person.
Collecting the data in a $n \times p$ matrix $X$. Each column corresponds to a particular variate (measurement). Each row corresponds to a sample (person).

### Sample mean [slide 91]

For each variate (measurement type) $j$:

$$ \bar{x}_j := 1/n \sum_{i=1}^n X_{ij} $$

where $n$ is the number of samples and $\bar{x}$ is the sample mean of each column of $X$. Using matrix notation:

$$ \bar{x} = 1/n X^T 1_n $$

where $1_n$ is column vector of ones with length $n$.
The $j$th entry in $\bar{x}$ is the sample mean of the $j$th column of the matrix $X$.

### Sample covariance matrix [slide 92]

The sample covariance matrix $S$ ($p \times p$) is given by:

$$ S_{jk} := 1/(n-1) \sum_{i=1}^n (X_{ij} - \bar{x}_j)(X_{ik} - \bar{x}_k) $$

In terms of the sample matrix $X$, the sample covariance matrix is given by:

$$ S = 1/(n-1) (X - 1_n \bar{x}^T)^T (X - 1_n \bar{x}^T) $$

### Sample mean and variance of the data [slide 93]

Recall that $d = X a$ is the change in the data (variates) for each sample in the direction of $a$. We want to find the direction $a$ in which the variance of the entries of $d$ is greatest.
Sample mean of $d = X a$:

$$ \bar{d} = 1/n 1_n^T X a = a^T \bar{x} $$

Sample variance of $d = X a$:

$$ q := 1/(n-1) \sum_{i=1}^n ((X a)_i - \bar{d})^2 = a^T S a $$

### Simple expression for the sample covariance [slide 94]

Adjusting $X$ to make it column centred, i.e. the columns of $X$ have zero mean, $\tilde{X} = X - 1_n \bar{x}^T$, gives a simple expression for the sample covariance matrix:

$$ (n-1)S = \tilde{X}^T \tilde{X} $$

and the sample variance of $d$ becomes

$$ q = 1/(n-1) a^T \tilde{X}^T \tilde{X} a $$

Introducing the centering matrix $C_n = I_{n \times n} - 1/n 1_n 1_n^T$, we have $\tilde{X} = C_n X$.

### Maximising variance [slide 95]

- Given the restriction $||a||_2 = 1$, we know that the sample variance $a^T S a$ is maximised when $a$ is aligned with the eigenvector associated with the largest eigenvalue of $S$.
- The eigenvalues of $(n-1)S$ are the singular values of $\tilde{X}$ (the multiplier $(n-1)$ is of no consequence since were interested in the relative sizes of the singular values).
- From $(n-1)S = \tilde{X}^T \tilde{X} = V \Sigma^T \Sigma V^T$, the matrix $V$ rotates the data to a 'new' coordinate system, ordered by decreasing variance.
- PCA is an eigen decomposition (diagonalisation) of $\tilde{X}^T \tilde{X}$.
- SVD algorithms are used to compute the singular values of $\tilde{X}$ and the eigenvectors $V$ without having to form $\tilde{X}^T \tilde{X}$.

See Jupyter notebooks for an application of PCA.

---

# Stochastic Processes

## Finite State-Space Markov Chains

### Overview

- This part of the course will focus of Markov Chains and related applications
- Course has three parts
  - finite-space Markov chains: 2 lectures
  - continuous state-space systems: 1.5 lectures
  - Monte Carlo Markov chains: 1.5 lectures
- Handouts available from 3M1 web-site
- Insufficient time to examine all topics in detail
  - useful starting point for 4th year modules

### Applications of Markov Chains

- Markov Chains under-pin the work in many areas
  - Google Page ranker (very large Markov chain!)
  - Information theory (Entropy)
  - Speech and Language Processing (acoustic models/language models)
  - Physics (statistical mechanics/thermodynamics)
  - Economics and finance (asset pricing)
  - Queueing theory

### Pushkin's