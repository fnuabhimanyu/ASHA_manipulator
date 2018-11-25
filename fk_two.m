function[Tf] = fk_two(links,a,alpha,d,theta)
syms F T
%L is the length of the links
%theta is the input angle
F=eye(4,4);
for i=1:1:links
T=[cos(theta(i)) -sin(theta(i))*cos(alpha(i)) sin(theta(i))*sin(alpha(i)) a(i)*cos(theta(i));
   sin(theta(i)) cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) a(i)*sin(theta(i));
    0            sin(alpha(i))                cos(alpha(i))               d(i);
    0            0                            0                           1];
F = F *T;
%x(i)=F(1,4);
%y(i)=F(2,4);
%z(i)=F(3,4);

end
Tf = F(1:3,1:3);
