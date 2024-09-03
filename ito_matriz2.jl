using Plots
using Distributions

function f(x)
    return r.*x.*(1 .-(a*x)./k) 
end

r = [1,1]
a = [1 1
    1 1  ]
b = [1, 1]*0.01
k = [60, 60]

tmax = 50
dt = 0.001
L = Int(tmax/dt)
tempo = collect(0:dt:tmax)
normal = Normal(0,sqrt(dt))
  
function ito_matriz(x1,x2)
    x = [x1,x2]
    xlist = Array{Float64}(undef,2, L+1)
    xlist[1,1], xlist[2,1] = x1,x2
    for t in 1:L
        random = [rand(normal),rand(normal)].*sqrt(dt)
        k1 = f(x)
        k2 = f(x .+ k1*dt/3)
        k3 = f(x .+ (-k1/3 + k2)*dt)
        k4 = f(x .+ (k1 - k2 + k3)*dt)
     
        dx = ((k1 + 3 .*k2 + 3 .*k3 + k4)./8 +0.5 .*x.*b.^2)*dt + b.*random .*x
        x = x + dx
        xlist[1,t+1], xlist[2,t+1] = x[1],x[2]
    end
    return xlist
end

result = stratmatriz(10,10)
x1,x2 = collect(result[1,:]),collect(result[2,:])
plot(tempo,x1)
p1 = plot!(tempo,x2,fmt = :png)

result = stratmatriz(10,10)
x1,x2 = collect(result[1,:]),collect(result[2,:])
plot(tempo,x1)
p2 = plot!(tempo,x2,fmt = :png)

l = @layout [a b ]
plot(p1, p2, layout = l, size=(1000,300))