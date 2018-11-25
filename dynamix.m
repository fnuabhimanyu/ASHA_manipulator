clc
clear
syms q1 q2 q3 q4 q5 qd1 qd2 qd3 qd4 qd5 qdd1 qdd2 qdd3 qdd4 qdd5 q6 qd6 qdd6 omg omgd vd vc m1 m2 m3 m4 m5 m6 F N f n tq
links = 6;
theta = [q1 q2 q3-sym(pi)/2 q4-sym(pi)/2 q5+sym(pi) 0];     % these 4 are column-indices of the respective JLP in A
d = [1 0 0 0 0 6];
a = [0 -2 -3 -4 5 0];
alpha = [sym(pi)/2 0 0 0 -sym(pi) 0];

% a = [0,0,0,0,5];
% alpha = [sym(pi)/2,sym(pi)/2,0,sym(pi)/2,sym(pi)/2];
% d = [0,0,5,5,0];
% theta = [sym(pi),(-sym(pi)/2)+q1,q2,sym(pi)+q3,(sym(pi)/2)+q4];

tp = [0,0,0,0,0,1];
qd = [qd1,qd2,qd3,qd4,qd5,qd6];
qdd = [qdd1,qdd2,qdd3,qdd4,qdd5,qdd6];
M = [m1 m2 m3 m4 m5 m6];
omg(1:3,1) = [0;0;0];
omgd(1:3,1) = [0;0;0];
vd(1:3,1) = [0;0;0];
f(1:3,7) = [0 0 0]';
n(1:3,7) = [0 0 0]';
r = [0 1 1.5 2 -2.5 0;
    -1/2 0 0 0 0 0;
    0 0 0 0 0 -3];

for i = 1:6
    
    T=[cos(theta(i)) -sin(theta(i))*cos(alpha(i)) sin(theta(i))*sin(alpha(i)) a(i)*cos(theta(i));
        sin(theta(i)) cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) a(i)*sin(theta(i));
        0            sin(alpha(i))                cos(alpha(i))               d(i);
        0            0                            0                           1];
    %P = ['T ', num2str(i)];
    
    dd = transpose(fk_two(i,a,alpha,d,theta))
    if tp(i) == 0
        omg(1:3,i+1) = T(1:3,1:3)'*(omg(1:3,i)+(transpose([0,0,qd(i)])));
        omgd(1:3,i+1) =  T(1:3,1:3)'*(omgd(1:3,i)+ (transpose([0,0,qdd(i)]))+cross(omg(1:3,i),(transpose([0,0,qdd(i)]))));
         vd(1:3,i+1) =  T(1:3,1:3)'*vd(1:3,i)+cross(omgd(1:3,i+1),(dd*T(1:3,4)))+ cross(omg(1:3,i+1),cross(omg(1:3,i+1),(dd*T(1:3,4))));
    else %tp[i] == 1
        omg(1:3,i+1) =  T(1:3,1:3)'*omg(1:3,i);
        omgd(1:3,i+1) =  T(1:3,1:3)'*omgd(1:3,i);
        vd(1:3,i+1) =  T(1:3,1:3)'*(vd(1:3,i) + transpose([0,0,qdd(i)])) + cross((2*omg(1:3,i+1)),( T(1:3,1:3)'*transpose([0,0,qd(i)]))) + cross(omgd(:,i+1),(dd*T(1:3,4)))+ cross(omg(:,i+1),cross(omg(:,i+1),(dd*T(1:3,4))));
    end
    vc(1:3,i) = vd(1:3,i+1) + cross(omg(1:3,i+1),(r(1:3,i))) +cross(omg(:,i+1),cross(omg(:,i+1),r(1:3,i)));
    
end
for j = links:-1:1
    
    if j == links
        Tf = eye(3);
    else
        Tf = fk(j+1,a,alpha,d,theta);
    end
    D = fk(j,a,alpha,d,theta);
    ddd = transpose(fk_two(j,a,alpha,d,theta));
    I = MoI(j);
    F(1:3,j) = M(j)*vc(j);
    N(1:3,j) = I*omgd(1:3,i+1)+ cross(omg(1:3,i+1),(I*omg(1:3,i+1)));
    f(1:3,j) = F(1:3,j) + Tf(1:3,1:3)*f(1:3,j+1);
    n(1:3,j) = Tf(1:3,1:3)*n(1:3,j+1) + cross((ddd*D(1:3,4)),Tf(1:3,1:3)*f(1:3,j+1)) + cross(((ddd*D(1:3,4)) + ddd*r(1:3,j)),F(1:3,j)) + N(1:3,j);
    if tp(j) == 0
        tq(j) = n(1:3,j)'*T(1:3,1:3)'*[0 0 1]'
    else
        tq(j) = f(1:3,j)'*T(1:3,1:3)'*[0 0 1]'
    end
end


