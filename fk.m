function [T] = fk(i,a,alpha,d,theta)

 T=[cos(theta(i)) -sin(theta(i))*cos(alpha(i)) sin(theta(i))*sin(alpha(i)) a(i)*cos(theta(i));
        sin(theta(i)) cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) a(i)*sin(theta(i));
        0            sin(alpha(i))                cos(alpha(i))               d(i);
        0            0                            0                           1];
end