#Problem 1 : unsteady 1D heat source conduction
#Method : Finite volume (explicit)
#Author : Chrislain Coubard
#Date : 27/05/2024

using PlotlyJS

function computeT(nt, nx, Tini)
    T = zeros(nt,nx)
    T[1,:] .= Tini
    for i = 1:(nt-1)
        for j = 1:nx
            if j == 1
                T[i+1,j] = Fo*T[i,j+1] + (1-Fo)*T[i,j]
            elseif j==nx
                T[i+1,j] = Fo*T[i, j-1] + (1-3*Fo)*T[i,j]+2*Fo*TL
            else
                T[i+1,j] = Fo*T[i,j+1] - (2*Fo-1)*T[i,j]+Fo*T[i,j-1]
            end
        end
    end
    return T
end

function Taugmented(T, init, TL, nt)
    TleftB = zeros(nt, 1)
    TleftB[1] = init
    TleftB[2:nt,:] .= T[2:nt,1]
    TrightB = zeros(nt,1)
    TrightB[1] = init
    TrightB[2:nt,:] .= TL
    Taugmented = [TleftB T TrightB]
    return Taugmented
end

rhocp = 1e7 #J/m^3/K
k=10 #W/m/K
a=k/rhocp #m^2/s

#Discretization
L = 0.02 #m
nx = 5 #Number of nodles (no boundaries)
dx = L/nx #Step Size

tf = 500  #s

#Minimum nt for stability

dt_min = dx^2/(2*a)
nt_min = 2*nx^2*a*tf/(L^2)+1
dt_min = tf / (nt_min-1)

#Number of time nodles including boundaries

nt = 101
dt = tf/(nt-1)

Fo = a*dt/dx^2

#Boundary conditions
TL = 0

#Initial conditions
Tini = 200

T = computeT(nt, nx, Tini)
Taug = Taugmented(T, Tini, TL, nt)

#Plotting the results

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
    z=Taug,x=t,y=x, colorscale = "Earth"),layout)

