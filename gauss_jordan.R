read_augmented_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  n <- as.integer(trimws(lines[1]))
  values <- scan(
    text = paste(lines[-1], collapse = " "),
    what = numeric(),
    quiet = TRUE
  )
  matrix(values, nrow = n, ncol = n + 1, byrow = TRUE)
}

gauss_jordan <- function(aug, tol = 1e-12) {
  n <- nrow(aug)

  for (k in seq_len(n)) {
    pivot_row <- which.max(abs(aug[k:n, k])) + k - 1

    if (abs(aug[pivot_row, k]) < tol) {
      stop("Matriz singular ou quase singular.")
    }

    if (pivot_row != k) {
      temp <- aug[k, ]
      aug[k, ] <- aug[pivot_row, ]
      aug[pivot_row, ] <- temp
    }

    aug[k, ] <- aug[k, ] / aug[k, k]

    for (i in seq_len(n)) {
      if (i == k) {
        next
      }

      factor <- aug[i, k]
      if (abs(factor) > tol) {
        aug[i, ] <- aug[i, ] - factor * aug[k, ]
      }
    }
  }

  list(x = aug[, n + 1], rref = aug)
}

residual_norm <- function(A, x, b) {
  r <- as.vector(b - A %*% x)
  list(vector = r, norm2 = sqrt(sum(r^2)))
}

refine_once <- function(A, b, x, tol = 1e-12) {
  res <- residual_norm(A, x, b)
  correction <- gauss_jordan(cbind(A, res$vector), tol = tol)$x
  x_new <- x + correction
  res_new <- residual_norm(A, x_new, b)

  list(
    x = x_new,
    correction = correction,
    before = res,
    after = res_new
  )
}

format_vector <- function(x, digits = 10) {
  paste(formatC(x, format = "fg", digits = digits), collapse = " ")
}

run_case <- function(path) {
  aug <- read_augmented_matrix(path)
  A <- aug[, -ncol(aug), drop = FALSE]
  b <- aug[, ncol(aug)]

  solution <- gauss_jordan(aug)$x
  refined <- refine_once(A, b, solution)

  before_norm <- refined$before$norm2
  after_norm <- refined$after$norm2
  ratio <- if (before_norm == 0) 1 else after_norm / before_norm

  cat("Arquivo:", path, "\n")
  cat("n =", nrow(A), "\n")
  cat("Solucao inicial:\n")
  cat(format_vector(solution), "\n")
  cat("Norma-2 do residuo antes do refinamento:", format(before_norm, digits = 16), "\n")
  cat("Solucao apos 1 iteracao de refinamento:\n")
  cat(format_vector(refined$x), "\n")
  cat("Norma-2 do residuo apos o refinamento:", format(after_norm, digits = 16), "\n")

  if (after_norm < before_norm) {
    cat("Analise: a norma do residuo diminuiu por um fator de", format(ratio, digits = 6), ".\n")
  } else if (after_norm > before_norm) {
    cat("Analise: a norma do residuo aumentou por um fator de", format(ratio, digits = 6), ".\n")
  } else {
    cat("Analise: a norma do residuo permaneceu inalterada.\n")
  }

  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)
files <- if (length(args) > 0) args else c("m1.in", "m2.in")

for (file in files) {
  run_case(file)
}
