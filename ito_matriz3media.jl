using Plots
using Distributions

r = [1,1,1]
a = [1 2 1
        1 1 2 
        2 1 1 ]
b = [0.1, 0.3, 0.03]
k = [60, 60, 60]
kinv = 1 ./k

tmax = 70
dt = 0.001
L = Int(tmax/dt)
tempo = collect(0:dt:tmax)
normal = Normal(0,sqrt(dt))
  N = 200 #n iteracoes
function stratmatriz(x1,x2,x3)
    
    xlist = zeros(Float64,3, L+1)
    xmlist = zeros(Float64,3, L+1)
    xvlist = zeros(Float64, 3, L+1)
    xlist[1,1], xlist[2,1], xlist[3,1], = x1,x2,x3
    for n in 1:N
        x = [x1,x2,x3]
        for t in 1:L
            random = [rand(normal),rand(normal),rand(normal)].*sqrt(dt)
    
            k1 = f(x)
            k2 = f(x .+ k1*dt/3)
            k3 = f(x .+ (-k1/3 + k2)*dt)
            k4 = f(x .+ (k1 - k2 + k3)*dt)
            
            dx = ((k1 + 3 .*k2 + 3 .*k3 + k4)./8 )*dt + x.*b.*random
            x = x + dx
            xlist[1,t+1], xlist[2,t+1], xlist[3,t+1] = xlist[1,t+1] + x[1],xlist[2,t+1]+x[2],xlist[3,t+1]+x[3]
            xvlist[1,t+1], xvlist[2,t+1], xvlist[3,t+1] = xvlist[1,t+1]+ x[1]^2,xvlist[2,t+1]+x[2]^2,xvlist[3,t+1]+x[3]^2
           
        end
    end
    xmlist = xmlist + xlist
    xvlist = (xvlist - (xmlist .* xmlist)/N)/N
    xmlist = xmlist / N
    xmlist[1,1], xmlist[2,1], xmlist[3,1], = x1,x2,x3
    return xmlist, xvlist
end


result, variance = stratmatriz(10,10,10) 

x1,x2,x4 = collect(result[1,:]),collect(result[2,:]),collect(result[3,:])
xs1,xs2,xs4 = collect(variance[1,:]),collect(variance[2,:]),collect(variance[3,:])
plot(tempo,x1)
plot!(tempo,x2)
#plot!(tempo,x3)
p1 = plot!(tempo,x4, fmt = :png, title = "média com ruído")  


b =0
result, variance= stratmatriz(10,10,10)
x1,x2,x4 = collect(result[1,:]),collect(result[2,:]),collect(result[3,:])
plot(tempo,x1)
plot!(tempo,x2)
#plot!(tempo,x3)
p2 = plot!(tempo,x4,fmt = :png, title = "média sem ruído")


l = @layout [a b ]
plot(p1, p2, layout = l, size=(1000, 300), fmt=:png, ylabel="Número de indivíduos", xlabel="Tempo", left_margin=10Plots.mm, bottom_margin=10Plots.mm)