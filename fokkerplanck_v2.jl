#teste
using Plots
using Distributions

function h(index, variable)
    return (r.*x[index].*(1 .-(a*x[index])./k))[variable]
end

    tmax = 6
    dt = 0.005

    x1max = 60
    x2max = 60
    
    dx1 = 0.4
    dx2 = 0.4
    
    x1 = collect(-dx1:dx1:x1max)
    x2 = collect(-dx2:dx2:x2max)
    
    L = Int(x1max/dx1)
    N = Int(tmax/dt)
    
    x = [[x1[i], x2[i]] for i = 1:L+2]
    
    xm = zeros(N)
    ym = zeros(N)
    
    espaco = [x1, x2] 
    r = [1,1]
    a = [1 1
        1 1  ]
    b = [1, 1]*0.01
    k = [60, 60]


function fokkerPlanck(x1i, x2i)
    
    P = zeros(L+2, L+2)
    G = zeros(L+2, L+2)
    px1 = Int(x1i/dx1)
    px2 = Int(x2i/dx2)
    P[px1,px2] = 1
    
     anim = @animate for t in 1:N
       for i in 2:L+1
            for j in 2:L+1
                G[i,j] = P[i,j] - dt*(h(i+1, 1)*P[i+1,j]-h(i-1, 1)*P[i-1,j])/2/dx1 +
                - dt*(h(j-1, 2)*P[i,j+1]-h(j-1, 2)*P[i,j-1])/2/dx2 +
                + dt*1/2/dx1^2*b[1]^2*(x[i+1][1]^2*P[i+1, j]+x[i-1][1]^2*P[i-1, j]-2*x[i][1]^2*P[i,j]) +
                + dt*1/2/dx2^2*b[2]^2*(x[j+1][2]^2*P[i, j+1]+x[j-1][2]^2*P[i, j-1]-2*x[j][2]^2*P[i,j])
                if G[i,j] < 0 
                    G[i, j] = 0 
                end
            end
        end
        P = copy(G)
        for i in 1:L+2
            for j in 1:L+2
                xm[t] = xm[t] + x1[i]*P[i,j]*dx1^2
                ym[t] = ym[t] + x2[j]*P[i,j]*dx2^2
                end
            end
        Plots.heatmap(espaco,espaco,P)
        end
        gif(anim, "fokkerplanck4.gif", fps=30)
        return P, xm, ym
    end

P, X, Y = fokkerPlanck(10, 10)

#sum(P)