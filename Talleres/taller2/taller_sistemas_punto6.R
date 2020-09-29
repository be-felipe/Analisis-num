A <- matrix(c(-8.1, -7, 6.123, -2, -1, 4, -3, -1, 0, -1, -5, 0.6, -1, 0.33, 6, 1/2), nrow=4, byrow=TRUE)

#Descomposicion LU
descLU <- lu( A )
lu_Desc <- expand(descLU)
L <- lu_Desc$L
U <- lu_Desc$U

print("Matriz L.")
print(L)
print("Matriz U.")
print(U)

print("Matriz Original A. ")
print(L %*% U)

#Descomposicion QR
qr_Desc <- qr( A ) 
Q <- qr.Q( qr_Desc )
R <- qr.R( qr_Desc )

print("Matriz Q.")
print(Q)
print("Matriz R.")
print(R)
print("Matriz Original A.")
print(Q %*% R)