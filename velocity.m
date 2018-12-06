clc
syms q1 q2 q3 q1d q2d
syms F R P Jw Jv J
D = [0 0 0 1]'
links = 2;
a = [ 2 2];
d = [ 0 0];
theta = [q1 q2 ];
alpha = [ 0 0];
Q = [q1 q2];
Qd = [q1d q2d];
for i = 1:1:links
    F(1:4,1:4,i+1) = fk_two(i,a,alpha,d,theta)
end
V1(1:4,1) = diff(F(1:4,1:4,2),Q(1))*[ 0 0 0 1]'*Qd(1);
W1(1:3,1) = [0 0 0]'+[1 0 0;0 1 0;0 0 1]*[0 0 1]'*Qd(1);
V1(1:4,2) = diff(F(1:4,1:4,3),Q(1))*[ 0 0 0 1]'*Qd(1)+diff(F(1:4,1:4,3),Q(2))*[ 0 0 0 1]'*Qd(2);
W1(1:3,2) = W1(1:3,1)+F(1:3,1:3,2)*[0 0 1]'*Qd(2);

V1
W1
    

    
