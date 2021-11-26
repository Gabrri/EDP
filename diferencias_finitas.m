%CODIGO GABRIELA CORTES MEJIA Y JUAN CAMILO RODRIGUEZ
clear all;
clc;

%Parametros iniciales
a=0;
b=1;
c=0;
d=2;
n=15;
m=15;

%Funcion
f=@(x,y) 2*x;

%Condiciones de frontera
uy_05=@(y) 2*y^2;

ux_05=@(x) x;

u0y=@(y) 0;
u0x=@(x) 0;

%Solucion Analitica
analitic = @(x,y) 0;
A_error = 0;
vector_sol_analitic = zeros((m-1)*(n-1),1);

h=(b-a)/n;
k=(d-c)/m;

vec_x = a+h:h:b-h;
vec_y = c+k:k:d-k;

w=zeros((m-1)*(n-1));
z=zeros((m-1)*(n-1),1);

lambda=h^2/k^2;
mu=2*(lambda+1);


for i=1:n-1
    for j=1:m-1
        xi= a+i*h;
        yj= c+j*k;
        values=i+(m-1-j)*(n-1);
        z(values)= z(values)+(-(h^2)*f(xi,yj));
        v_1= i+(m-1-j)*(n-1);
        w(values,v_1)=mu;

        if A_error == 1
            vector_sol_analitic(values) = analitic(xi,yj);
        end

        if i+1 ~= n
            v_2=(i+1)+(m-1-j)*(n-1);
            w(values,v_2)=-1;
        else
            z(values)=z(values)+uy_05(yj);
        end

        if i-1~=0
            v_3= i-1 +(m-1-j)*(n-1);
            w(values,v_3)= -1;
        else
            z(values)=z(values)+u0y(yj);
        end

        if j+1~=m
            v_4=i+(m-1-(j+1))*(n-1);
            w(values,v_4)=-lambda;
        else
            z(values)=z(values)+lambda*ux_05(xi);
        end

        if j-1~=0
            v_5 = i+(m-1-(j-1))*(n-1);
            w(values,v_5) = -lambda;
        else
            z(values) = z(values)+lambda*u0x(xi);
        end

    end
end

%Metodo del gradiente conjugado
sol = bicg(w,z);

vec_error = vector_sol_analitic - sol;

M_sol = zeros((m-1),(n-1));
M_error = zeros((m-1),(n-1));
for i=1:(m-1)
        for j=1:(n-1)
            val = i + (n-1-j)*(m-1);
            M_sol(i,j) = sol(val);
            M_error(i,j) = vec_error(val);
        end
end

[x,y] = meshgrid(vec_x,vec_y);

surf(x,y,M_sol)
colorbar
title("Solucion con Diferencias Finitas")
