function gamma1(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
    return (r1*x1^2)/k1 + r1*x1*a12*x2/k1 - r1*x1 + 2*x1*b1^2
end

function gamma2(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
    return (r2*x2^2)/k2 + r2*x2*a21*x1/k2 - r2*x2 + 2*x2*b2^2
end

function delta1(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
    return (b1^2 * x1^2)/2
end    

function delta2(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
    return (b2^2 * x2^2)/2
end

function omega(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
    return 2*r1*x1/k1 + 2*r2*x2/k2 + r1*x2*a12/k1 + r2*x1*a21/k2 - r1 - r2 + b1^2 + b2^2
end

xmax = 60
tmax = 6
dx = 0.4
dt = 0.005

N = Int(tmax/dt)
L = Int(xmax/dx)

tempo = collect(0:dt:tmax)
espaco = collect(0:dx:xmax)

k1, k2, r1, r2, a12, a21, b1, b2 = 60, 60, 1, 1, 1, 1, 0.01, 0.01

function fokkerplank(x10,x20)
    P = zeros(L+1, L+1)
    G = zeros(L+1, L+1)
    px1 = Int(x10/dx)
    px2 = Int(x20/dx)
    P[px1,px2] = 1
    
    @gif for t in 1:N
       for i in 2:L
            x1 = espaco[i]
            for j in 2:L
                x2 = espaco[j]
                
                g1 = gamma1(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
                g2 = gamma2(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
                d1 = delta1(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
                d2 = delta2(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
                w = omega(x1,x2,r1,r2,k1,k2,a12,a21,b1,b2)
                G[i,j] = P[i,j]*(1 - 2*d1*dt/dx^2 - 2*d2*dt/dx^2 + w*dt) + P[i+1,j]*(d1*dt/dx^2 + 0.5*g1*dt/dx) + P[i-1,j]*(d1*dt/dx^2 - 0.5*g1*dt/dx) + P[i,j+1]*(d2*dt/dx^2 + 0.5*g2*dt/dx) + P[i,j-1]*(d2*dt/dx^2 - 0.5*g2*dt/dx)
                if G[i,j] < 0 
                    G[i, j] = 0 
                end
            end
        end
        P = copy(G)
        Plots.heatmap(espaco,espaco,P)
    end 
    return P
end

p = fokkerplank(10,10)