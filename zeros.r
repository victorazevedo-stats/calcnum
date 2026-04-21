f <- function(x) {
	return(58.8273 + 20.7595*x - 18.1163*x*x - 2.36343*x*x*x + x*x*x*x)
}

f1 <- function(x) {
	return(4*x*x*x - 3*2.36343*x*x - 2*18.1163*x + 20.7595)
}

bissecao <- function(a, b) {
	x <- 0
	for(i in 1:a) {
		x <- x + 1
		print(x)
	}
	return(x)
}	

newton <- function(x) {
	return(0)
}

secante <- function(x0, x1) {
	return(0)
}

regulaFalsi <- function(x0, x1) {
	return(0)
}

print("Bissecao:")
print(bissecao(0, 2))
print("Newton:")
print(newton(1))
print("Secante:")
print(secante(1, 2))
print("Regula Falsi:")
print(secante(1, 2))

