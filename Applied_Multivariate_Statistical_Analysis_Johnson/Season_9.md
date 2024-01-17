# Excercises of Season 9

#### Ex. iv.

9.4. Given $\rho$ and $\psi$ in Exercise 9.1 and an m = 1 factor model,
calculate the reduced correlation matrix $\tilde{\rho} = \rho - \psi$
and the principal factor solution for the loading matrix L. Is the
result consistent with the information in Exercise 9.1? Should it be?

###### Soloution

``` r
rho <- matrix(c(1, 0.63, 0.45, 0.63, 1, 
                .35, .45, .35, 1), 3, 3, byrow = T)
rho
```

         [,1] [,2] [,3]
    [1,] 1.00 0.63 0.45
    [2,] 0.63 1.00 0.35
    [3,] 0.45 0.35 1.00

``` r
psi <- diag(c(0.19, 0.51, 0.75), 3)
psi
```

         [,1] [,2] [,3]
    [1,] 0.19 0.00 0.00
    [2,] 0.00 0.51 0.00
    [3,] 0.00 0.00 0.75

``` r
L <- matrix(c(0.9, 0.7, 0.5), 3, 1)
ll <- L %*% t(L)

le <- eigen(rho)$values
ve <- eigen(rho)$vectors
ress <- rho - psi

all(ress - ll < 1e-10)
```

    [1] TRUE

<br>

------------------------------------------------------------------------

------------------------------------------------------------------------

### EX. 9.8

9.8. (Unique but improper solution: Heywood case.) Consider an m = 1
factor model for the population with covariance matrix

$$
\left[\begin{array}{ccc}
1 & 0.4 & 0.9 \\
0.4 & 1 & 0.7 \\
0.9 & 0.7 & 1
\end{array}\right]
$$

Show that there is a unique choice of $L$ and $\Psi$ with
$\Sigma = L \times L^T + \Psi$, but that $\psi < 0$, so the choice is
not admissible

###### Soloution:

$$
\begin{aligned}
& \Sigma = LL' + \Psi \quad \text{for} ~m = 1 \implies \\
& \Sigma = \left[\begin{array}{ccc} 
1 = L_{11}^2 + \psi_1 & 0.4 = L_{11} \times L_{21} & 0.9 = L_{11} \times L_{31} \\
0.4 = L_{11} \times L_{21} &  1 = L_{21}^2 + \psi_2 & 0.7 = L_{21}\times L_{31} \\
0.9 = L_{11} \times L_{31} & 0.7 = L_{21}\times L_{31} & 1 = L_{31}^2 + \psi_3
\end{array}\right] \implies \\
& \text{We Have:}\quad \frac{L_{11}}{L_{21}} = \frac{0.9}{0.7} ~\text{and}~L_{11} \times L_{21} = 0.4, ~ \text{So} \quad \\
& L_{11}^2 = \frac{0.9}{0.7} \times 0.4 ~ \text{and} ~ L_{11} = \pm 0.717 \implies\\
& \text{Thus}\quad L_{21} = \pm 0.558, \quad \text{Finally}\\
& \text{from}\quad 0.9 = L_{11}\times L_{31} \implies \\
& \text{We Have:}\quad L_{31} = \pm \frac{0.9}{0.717} = \pm 1.255
\end{aligned}
$$

Note all the loadings must be of the same sign because all the
covariances are positive. We have

$$
\begin{aligned}
& L\times L^T = \left[\begin{array}{c} 0.717 \\ 0.558 \\ 1.255\end{array}\right] \times \left[\begin{array}{r}0.717 & 0.558 & 1.255\end{array}\right] = \\
& \left[\begin{array}{ccc} 0.514 & 0.4 & 0.9 \\ 0.4 & 0.3111 & 0.7 \\
0.9 & 0.7 & 1.575 \end{array}\right] \implies \\
& \text{So} \quad \psi_3 = 1  - 1.575 = -0.575 < 0 \text{which is inadmissible as a variance.}
\end{aligned}
$$