#Problem 1 : unsteady 1D heat source conduction
#Method : Finite volume (implicit)
#Author : Chrislain Coubard
#Date : 27/05/2024

using PlotlyJS
using LinearAlgebra

function maketridiag(a1, aEnd, aN, b,n)
    MainDiag = fill(aN,n)
    MainDiag[1] = a1
    MainDiag[n] = aEnd
    OtherDiags = fill(b,n-1)
    Mat = Tridiagonal(OtherDiags, MainDiag, OtherDiags)
    return Mat
end

function SolveT(nx, nt, A, B, Tini)
    T = zeros(nt,nx)
    T[1,:] .= Tini
    for i = 1:(nt-1)
        T[i+1,:] = A\B
        B[1:(nx-1)] = T[i+1,1:(nx-1)]
        B[nx] = T[i+1,nx]+2*Fo*TL
    end
    return T
end

function Taugmented(T, Tini, TL,nt)
    Tleft = zeros(nt)
    Tleft[1] = Tini
    Tleft[2:end] = T[2:end,1]
    Tright = zeros(nt)
    Tright[1] = Tini
    Tright[2:end] .= TL
    Taugmented = [Tleft T Tright]
    return Taugmented
end

t1 = time()
#Problem : Backward Euler 1st order difference scheme
rhocp = 1e7 #J/m3/K
k = 10 #W/m/K
a = k/rhocp #m2/s

#Discretization

L = 0.02 #m
nx = 5 #Number of x nodles (no boundarie included)
dx = L/nx #Step size
 
tf = 500 #s

#Number of time nodles including boundaries
nt = 101 
dt = tf/(nt-1) #Step size

#Fourier Number
Fo = a*dt/dx^2

# A = zeros(nx,nx)

#Boundary conditions
TL = 0 #°C

#Initial condition
Tini = 200 #°C
B = zeros(nx,1)
B[1:nx-1] .= Tini
B[nx] = Tini+2*Fo*TL

T[1,:] .= Tini

#remplissage matrice A

# for i = 1:nx
#     if i==1
#         A[i,i] = 1 + Fo
#         A[i,i+1] = -Fo
#     elseif i==nx
#         A[i,i] = 1+3*Fo
#         A[i,i-1] = -Fo
#     else
#         A[i,i] = 1+2*Fo
#         A[i, i-1] = -Fo
#         A[i, i+1] = -Fo
#     end
# end

#Make tridiagMatrix A
a1 = 1 + Fo
aEnd = 1 + 3*Fo
aN = 1 + 2*Fo
b = - Fo
A = maketridiag(a1,aEnd,aN,b,nx)

#Fill T matrix and add boudaries conditions
T = SolveT(nx, nt, A, B, Tini)
T = Taugmented(T, Tini, TL, nt)

t2 = time() - t1
println("Elapsed time ", t2, " seconds")

#Plot using PlotlyJS for interactivity
x = dx/2:dx:L-dx/2
x = vcat(0,x,L)

t = 0 : dt : tf

layout = Layout(
    title = "1D heat diffusion",
    scene = attr(
    xaxis_title = "Temps (s)",
    yaxis_title = "Position (m)",
    zaxis_title = "Température (°C)")
    )


plot(surface(
    contours = attr(
        x=attr(show=true,start=0,size = dt*10),
        x_end=tf,
        y=attr(show=true,start=dx/2,size = dx),
        y_end =L), 
    z=T,x=t,y=x, colorscale = "Earth"),layout)
