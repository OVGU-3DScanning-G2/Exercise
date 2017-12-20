function [Center, Radius] = sphere_fit(x,y,z)
[N,dum]=size(x);

Sx = sum(x); 

Sy = sum(y);

Sz = sum(z);


Sxx = sum(x.*x);
Sxy = sum(x.*y);

Syy = sum(y.*y);
Syz = sum(y.*z);

Szz = sum(z.*z);
Sxz = sum(x.*z);


Sxxx = sum(x.*x.*x);
Sxxy = sum(x.*x.*y);
Sxyy = sum(x.*y.*y);

Syyy = sum(y.*y.*y);
Syyz = sum(y.*y.*z);
Syzz = sum(y.*z.*z);

Szzz = sum(z.*z.*z);
Sxzz = sum(x.*z.*z);
Sxxz = sum(x.*x.*z);


A1 = Sxx +Syy +Szz;

a=2*Sx*Sx-2*N*Sxx;
b=2*Sx*Sy-2*N*Sxy;
c=2*Sx*Sz-2*N*Sxz;
d= -N * (Sxxx +Sxyy +Sxzz) + A1 * Sx;

e = 2*Sx*Sy-2*N*Sxy;
f = 2*Sy*Sy-2*N*Syy;
g = 2*Sy*Sz-2*N*Syz;
h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;

j = 2*Sx*Sz-2*N*Sxz;
k = 2*Sy*Sz-2*N*Syz;
l = 2*Sz*Sz-2*N*Szz;
m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;

delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);
Center = zeros(3,1);
Center(1) = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
Center(2) = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
Center(3) = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;

Radius = sqrt(Center(1)^2+Center(2)^2+Center(3)^2+(A1-2*(Center(1)*Sx+Center(2)*Sy+Center(3)*Sz))/N);
end