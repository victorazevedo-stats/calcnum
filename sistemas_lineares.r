############################################################
# CÁLCULO NUMÉRICO - SISTEMAS LINEARES
# Códigos em R baseados nos tópicos dos slides
#
# Tópicos:
# 1) Representação matricial de sistemas lineares
# 2) Classificação de sistemas
# 3) Sistemas triangulares e substituição retroativa
# 4) Eliminação de Gauss
# 5) Eliminação de Gauss com pivotação parcial
# 6) Refinamento iterativo
# 7) Normas de vetores
# 8) Normas de matrizes
# 9) Método de Jacobi
# 10) Método de Gauss-Seidel
# 11) Critério de convergência (diagonal dominante)
# 12) Critérios de parada
# 13) Número de condicionamento
# 14) Decomposição LU
# 15) Determinante via LU
# 16) Matriz inversa
############################################################


############################################################
# 1) REPRESENTAÇÃO MATRICIAL DE UM SISTEMA LINEAR
# Ax = b
############################################################

A <- matrix(c(2, -1,  2,
             -2, -4, -1,
              2,  3,  1), 
            nrow = 3, byrow = TRUE)

b <- c(3, 5, 2)

cat("Matriz A:\n")
print(A)

cat("\nVetor b:\n")
print(b)

cat("\nMatriz aumentada [A|b]:\n")
print(cbind(A, b))


############################################################
# 2) CLASSIFICAÇÃO DE SISTEMAS
# Possível determinado, possível indeterminado, impossível
# usando postos (rank)
############################################################

classificar_sistema <- function(A, b, tol = 1e-10) {
  Ab <- cbind(A, b)
  posto_A  <- qr(A, tol = tol)$rank
  posto_Ab <- qr(Ab, tol = tol)$rank
  n <- ncol(A)
  
  if (posto_A == posto_Ab && posto_A == n) {
    return("Sistema possível e determinado (solução única)")
  } else if (posto_A == posto_Ab && posto_A < n) {
    return("Sistema possível e indeterminado (infinitas soluções)")
  } else {
    return("Sistema impossível (sem solução)")
  }
}

cat("\nClassificação do sistema A x = b:\n")
print(classificar_sistema(A, b))


############################################################
# 3) SISTEMAS TRIANGULARES
# Resolução por substituição retroativa
# Os slides destacam que sistemas triangulares são resolvidos
# por substituições retroativas.
############################################################

# Função para sistema triangular superior Ux = b
substituicao_retroativa <- function(U, b) {
  n <- nrow(U)
  x <- numeric(n)
  
  for (i in n:1) {
    if (abs(U[i, i]) < 1e-14) {
      stop("Pivô nulo encontrado na substituição retroativa.")
    }
    
    if (i == n) {
      soma <- 0
    } else {
      soma <- sum(U[i, (i+1):n] * x[(i+1):n])
    }
    
    x[i] <- (b[i] - soma) / U[i, i]
  }
  
  return(x)
}

# Exemplo dos slides
U <- matrix(c(3, -2,  3,
              0,  2, -1,
              0,  0,  2),
            nrow = 3, byrow = TRUE)

bU <- c(5, 3, 2)

cat("\nSistema triangular superior:\n")
print(U)
cat("\nSolução por substituição retroativa:\n")
print(substituicao_retroativa(U, bU))


############################################################
# 4) ELIMINAÇÃO DE GAUSS
# Transforma A em matriz triangular superior
############################################################

gauss_sem_pivotamento <- function(A, b) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  
  Ab <- cbind(A, b)
  
  for (k in 1:(n-1)) {
    if (abs(Ab[k, k]) < 1e-14) {
      stop("Pivô nulo encontrado. Use pivotação parcial.")
    }
    
    for (i in (k+1):n) {
      m <- Ab[i, k] / Ab[k, k]   # multiplicador
      Ab[i, k:(n+1)] <- Ab[i, k:(n+1)] - m * Ab[k, k:(n+1)]
    }
  }
  
  U <- Ab[, 1:n]
  c_final <- Ab[, n+1]
  x <- substituicao_retroativa(U, c_final)
  
  list(U = U, b_mod = c_final, x = x, Ab_final = Ab)
}

cat("\nEliminação de Gauss sem pivotamento:\n")
resultado_gauss <- gauss_sem_pivotamento(A, b)
print(resultado_gauss$Ab_final)
cat("\nSolução:\n")
print(resultado_gauss$x)


############################################################
# 5) ELIMINAÇÃO DE GAUSS COM PIVOTAÇÃO PARCIAL
# Nos slides, a pivotação parcial escolhe o maior elemento
# em módulo na coluna atual, da diagonal para baixo.
############################################################

gauss_pivotacao_parcial <- function(A, b) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  
  Ab <- cbind(A, b)
  
  for (k in 1:(n-1)) {
    # Escolha do pivô: maior valor absoluto na coluna k
    linha_pivo <- which.max(abs(Ab[k:n, k])) + k - 1
    
    if (abs(Ab[linha_pivo, k]) < 1e-14) {
      stop("Sistema singular ou quase singular.")
    }
    
    # Troca de linhas, se necessário
    if (linha_pivo != k) {
      temp <- Ab[k, ]
      Ab[k, ] <- Ab[linha_pivo, ]
      Ab[linha_pivo, ] <- temp
    }
    
    # Eliminação
    for (i in (k+1):n) {
      m <- Ab[i, k] / Ab[k, k]
      Ab[i, k:(n+1)] <- Ab[i, k:(n+1)] - m * Ab[k, k:(n+1)]
    }
  }
  
  U <- Ab[, 1:n]
  c_final <- Ab[, n+1]
  x <- substituicao_retroativa(U, c_final)
  
  list(U = U, b_mod = c_final, x = x, Ab_final = Ab)
}

# Exemplo inspirado nos slides de pivotação parcial
A_piv <- matrix(c( 0.12, -0.16,  0.08, -0.64,
                   0.35,  0.74,  0.57, -0.54,
                  -0.92,  0.07,  0.96,  0.40,
                  -0.15, -0.29, -0.74,  0.99),
                nrow = 4, byrow = TRUE)

b_piv <- c(0.31, 0.30, 0.18, -0.93)

cat("\nEliminação de Gauss com pivotação parcial:\n")
resultado_piv <- gauss_pivotacao_parcial(A_piv, b_piv)
print(resultado_piv$Ab_final)
cat("\nSolução:\n")
print(resultado_piv$x)


############################################################
# 6) REFINAMENTO ITERATIVO
# Slides: r = b - A x_bar
# Resolve A y = r
# Nova solução: x_ref = x_bar + y
############################################################

refinamento_iterativo <- function(A, b, x_aprox, n_iter = 3) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  x <- as.numeric(x_aprox)
  
  for (k in 1:n_iter) {
    r <- b - A %*% x
    y <- solve(A, r)      # resolve A y = r
    x <- x + as.vector(y)
    
    cat("\nIteração de refinamento", k, "\n")
    cat("Resíduo:\n")
    print(as.vector(r))
    cat("Correção y:\n")
    print(as.vector(y))
    cat("Nova aproximação:\n")
    print(x)
  }
  
  return(x)
}

# Aproximação inicial
x_aprox <- c(1, 1, 1)

cat("\nRefinamento iterativo:\n")
x_refinado <- refinamento_iterativo(A, b, x_aprox, n_iter = 2)


############################################################
# 7) NORMAS DE VETORES
# Slides: norma 1, norma 2 e norma infinito
############################################################

norma_1_vetor <- function(v) {
  sum(abs(v))
}

norma_2_vetor <- function(v) {
  sqrt(sum(v^2))
}

norma_inf_vetor <- function(v) {
  max(abs(v))
}

v <- c(1, 3, -2)

cat("\nNormas do vetor v = (1, 3, -2):\n")
cat("Norma 1:", norma_1_vetor(v), "\n")
cat("Norma 2:", norma_2_vetor(v), "\n")
cat("Norma infinito:", norma_inf_vetor(v), "\n")


############################################################
# 8) NORMAS DE MATRIZES
# Slides: ||A||_1, ||A||_2, ||A||_inf
# Em R, norm(A, type = "O"), "F", "I", "2"
# Para a norma 1 usamos type = "O" (one norm)
############################################################

norma_1_matriz <- function(A) {
  max(colSums(abs(A)))
}

norma_inf_matriz <- function(A) {
  max(rowSums(abs(A)))
}

norma_2_matriz <- function(A) {
  sv <- svd(A)$d
  max(sv)
}

A_ex <- matrix(c(1, 3, -2, 3, 2,
                 2, -2, -4, 0, 1,
                 1, -2, -5, 2, 2,
                 5, 2, -8, 2, 0,
                 0, 0, 1, -1, 0),
               nrow = 5, byrow = TRUE)

cat("\nNormas da matriz A_ex:\n")
cat("Norma 1:", norma_1_matriz(A_ex), "\n")
cat("Norma infinito:", norma_inf_matriz(A_ex), "\n")
cat("Norma 2:", norma_2_matriz(A_ex), "\n")


############################################################
# 9) MÉTODO DE JACOBI
# Slides: usa somente os valores da iteração anterior
############################################################

jacobi <- function(A, b, x0 = NULL, tol = 1e-8, max_iter = 100) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  
  if (is.null(x0)) x0 <- rep(0, n)
  
  x_old <- x0
  x_new <- numeric(n)
  
  for (k in 1:max_iter) {
    for (i in 1:n) {
      soma <- sum(A[i, -i] * x_old[-i])
      x_new[i] <- (b[i] - soma) / A[i, i]
    }
    
    erro_rel <- norm(x_new - x_old, type = "2") / max(norm(x_new, type = "2"), 1e-15)
    
    cat("\nJacobi - iteração", k, "\n")
    print(x_new)
    cat("Erro relativo:", erro_rel, "\n")
    
    if (erro_rel <= tol) {
      return(list(x = x_new, iter = k, convergiu = TRUE))
    }
    
    x_old <- x_new
  }
  
  return(list(x = x_new, iter = max_iter, convergiu = FALSE))
}

A_j <- matrix(c(10, 3, -2,
                2, 8, -1,
                1, 1, 5),
              nrow = 3, byrow = TRUE)

b_j <- c(57, 20, -4)

cat("\nMétodo de Jacobi:\n")
res_jacobi <- jacobi(A_j, b_j, x0 = c(0, 0, 0), tol = 1e-6, max_iter = 50)
cat("\nResultado final de Jacobi:\n")
print(res_jacobi)


############################################################
# 10) MÉTODO DE GAUSS-SEIDEL
# Slides: usa imediatamente os valores atualizados na mesma iteração
############################################################

gauss_seidel <- function(A, b, x0 = NULL, tol = 1e-8, max_iter = 100) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  
  if (is.null(x0)) x0 <- rep(0, n)
  
  x <- x0
  
  for (k in 1:max_iter) {
    x_old <- x
    
    for (i in 1:n) {
      soma1 <- if (i > 1) sum(A[i, 1:(i-1)] * x[1:(i-1)]) else 0
      soma2 <- if (i < n) sum(A[i, (i+1):n] * x_old[(i+1):n]) else 0
      x[i] <- (b[i] - soma1 - soma2) / A[i, i]
    }
    
    erro_rel <- norm(x - x_old, type = "2") / max(norm(x, type = "2"), 1e-15)
    
    cat("\nGauss-Seidel - iteração", k, "\n")
    print(x)
    cat("Erro relativo:", erro_rel, "\n")
    
    if (erro_rel <= tol) {
      return(list(x = x, iter = k, convergiu = TRUE))
    }
  }
  
  return(list(x = x, iter = max_iter, convergiu = FALSE))
}

cat("\nMétodo de Gauss-Seidel:\n")
res_gs <- gauss_seidel(A_j, b_j, x0 = c(0, 0, 0), tol = 1e-6, max_iter = 50)
cat("\nResultado final de Gauss-Seidel:\n")
print(res_gs)


############################################################
# 11) CRITÉRIO SUFICIENTE DE CONVERGÊNCIA
# Diagonal estritamente dominante
# Slides: vale tanto para Jacobi quanto para Gauss-Seidel
############################################################

eh_diagonal_dominante <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  for (i in 1:n) {
    if (abs(A[i, i]) <= sum(abs(A[i, -i]))) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

cat("\nA matriz A_j é diagonal estritamente dominante?\n")
print(eh_diagonal_dominante(A_j))


############################################################
# 12) CRITÉRIOS DE PARADA
# Slides:
# - por número máximo de iterações
# - por erro relativo
############################################################

criterio_parada_erro <- function(xk, xk1, eps = 1e-6) {
  erro <- norm(xk - xk1, type = "2") / max(norm(xk, type = "2"), 1e-15)
  return(erro <= eps)
}

xk  <- c(1.0001, 2.0001)
xk1 <- c(1.0000, 2.0000)

cat("\nCritério de parada por erro satisfeito?\n")
print(criterio_parada_erro(xk, xk1, eps = 1e-3))


############################################################
# 13) NÚMERO DE CONDICIONAMENTO
# Slides: cond(A) = ||A|| * ||A^{-1}||
############################################################

A_cond <- matrix(c(1, 0.99,
                   0.99, 0.98),
                 nrow = 2, byrow = TRUE)

cat("\nNúmero de condicionamento da matriz A_cond (norma 2):\n")
print(kappa(A_cond))

# Outro exemplo dos slides
A_cond2 <- matrix(c(1, 0.99,
                    0.99, 0.48),
                  nrow = 2, byrow = TRUE)

cat("\nNúmero de condicionamento da matriz A_cond2 (norma 2):\n")
print(kappa(A_cond2))


############################################################
# 14) DECOMPOSIÇÃO LU
# Slides: A = LU, com L triangular inferior e U superior
############################################################

decomposicao_LU <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  L <- diag(n)
  U <- matrix(0, n, n)
  
  for (i in 1:n) {
    # Calcula a linha i de U
    for (k in i:n) {
      U[i, k] <- A[i, k] - sum(L[i, 1:(i-1)] * U[1:(i-1), k])
    }
    
    # Calcula a coluna i de L
    for (k in (i+1):n) {
      if (abs(U[i, i]) < 1e-14) {
        stop("Pivô nulo encontrado na decomposição LU.")
      }
      L[k, i] <- (A[k, i] - sum(L[k, 1:(i-1)] * U[1:(i-1), i])) / U[i, i]
    }
  }
  
  list(L = L, U = U)
}

A_lu <- matrix(c(5, 2, 1,
                 3, 1, 4,
                 1, 1, 3),
               nrow = 3, byrow = TRUE)

res_lu <- decomposicao_LU(A_lu)

cat("\nMatriz L:\n")
print(res_lu$L)

cat("\nMatriz U:\n")
print(res_lu$U)

cat("\nProduto L %*% U:\n")
print(res_lu$L %*% res_lu$U)


############################################################
# 15) DETERMINANTE VIA LU
# Slides: det(A) = det(L) * det(U)
# Se diagonal de L é 1, então det(L) = 1
# Logo, det(A) = produto da diagonal de U
############################################################

det_via_lu <- function(A) {
  lu <- decomposicao_LU(A)
  prod(diag(lu$L)) * prod(diag(lu$U))
}

cat("\nDeterminante de A_lu via LU:\n")
print(det_via_lu(A_lu))

cat("\nDeterminante de A_lu via função det():\n")
print(det(A_lu))


############################################################
# 16) MATRIZ INVERSA
# Slides: A * A^{-1} = I
# Em R, pode-se usar solve(A)
############################################################

cat("\nMatriz inversa de A_lu:\n")
A_inv <- solve(A_lu)
print(A_inv)

cat("\nVerificação A %*% A^{-1}:\n")
print(A_lu %*% A_inv)


############################################################
# 17) RESOLVENDO SISTEMA COM solve()
# Apenas como comparação com os métodos implementados
############################################################

cat("\nSolução exata/computacional de A x = b com solve():\n")
print(solve(A, b))


############################################################
# 18) EXEMPLO EXTRA: RESOLVER LU EM DUAS ETAPAS
# Ly = b e depois Ux = y
############################################################

substituicao_progressiva <- function(L, b) {
  n <- nrow(L)
  y <- numeric(n)
  
  for (i in 1:n) {
    if (i == 1) {
      soma <- 0
    } else {
      soma <- sum(L[i, 1:(i-1)] * y[1:(i-1)])
    }
    y[i] <- (b[i] - soma) / L[i, i]
  }
  
  return(y)
}

resolver_por_lu <- function(A, b) {
  lu <- decomposicao_LU(A)
  y <- substituicao_progressiva(lu$L, b)
  x <- substituicao_retroativa(lu$U, y)
  
  list(L = lu$L, U = lu$U, y = y, x = x)
}

cat("\nResolvendo A_lu x = b via LU:\n")
b_lu <- c(10, 13, 7)
res_sis_lu <- resolver_por_lu(A_lu, b_lu)

cat("Vetor intermediário y:\n")
print(res_sis_lu$y)

cat("Solução x:\n")
print(res_sis_lu$x)