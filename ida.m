%% Inverse dynamic analysis of slider-crank mechanism
function ida()
m1 = 1;m2 = 2;m3 = 1;g = 9.81;
l1 = 0.4;l2 = 0.5;l3 =0.4;
J1 = m1*l1^2/12;J2 = m2*l2^2/12;J3 = m3*l3^2/12;
M = diag([m1 m1 J1 m2 m2 J2 m3 m3 J3]);
h = [0 -m1*g 0 0 -m2*g 0 0 -m3*g 0]';
load sckadata
N = size(t,2);
torque = zeros(1,N);
for i = 1:N
   [Phi,D] = constraints(t(i),pcoordsall(:,i));
   D_driver = [0 0 1 0 0 0 0 0 0];
   D_all = [D' D_driver'];
   rhs = M*acoordsall(:,i)-h;
   lambda = D_all\rhs;
   h_driver = D_driver'*lambda(9);% Direction * Magnitude
   h_driver';
   torque(i) = h_driver(3); % Finding Torque required at crank
end
plot(t,torque)
title('Input Torque');
xlabel('time');
ylabel('Torque');
end

%% Finding phi and system jacobian 

function [Phi,D] = constraints(t,q)
l1 = 0.4;l2 = 0.5;l3 = 0.4;
x1 = q(1);y1 = q(2);phi1 = q(3);
x2 = q(4);y2 = q(5);phi2 = q(6);
x3 = q(7);y3 = q(8);phi3 = q(9);

r1 = [x1 y1]';r2 = [x2 y2]';r3 = [x3 y3]';% w.r.to global frame
s_1_O = [-l1/2 0]';s_1_A = [l1/2 0]';% (Small letter's') Local position vectors for body 1 at joints  o and A w.r.to from body coordinate frame
s_2_A = [-l2/2 0]';s_2_B = [l2/2 0]';
s_3_B = [-l3/2 0]';s_3_C = [l3/2 0]';   

S_1_O = A(phi1)*s_1_O;S_1_A = A(phi1)*s_1_A; % (capital 'S') Global position vectors
S_2_A = A(phi2)*s_2_A;S_2_B = A(phi2)*s_2_B;
S_3_B = A(phi3)*s_3_B;S_3_C=A(phi3)*s_3_C;

Phi = [r1+S_1_O;  % phi is positin constraints from each joint
       r2+S_2_A-r1-S_1_A;
       r3+S_3_B-r2-S_2_B;
      [0.05 0]'-r3-S_3_C  
       ];
S_1_O_r = [0 -1;1 0]*S_1_O;S_1_A_r = [0 -1;1 0]*S_1_A;% Rotation of global vestors i.e, S rotated by 90 degrees of joint 'O' on body '1' = A*S_1_o 
S_2_A_r = [0 -1;1 0]*S_2_A;S_2_B_r = [0 -1;1 0]*S_2_B;
S_3_B_r = [0 -1;1 0]*S_3_B;S_3_C_r = [0 -1;1 0]*S_3_C;

D = [eye(2) S_1_O_r zeros(2,6);  % system jacobian matrix size (nc*nv)=(8*9)
    -eye(2) -S_1_A_r eye(2) S_2_A_r zeros(2,3);
    zeros(2,3) -eye(2) -S_2_B_r eye(2) S_3_B_r ;
    zeros(2,6) -eye(2) -S_3_C_r ;
    ];
end

function output = A(phi)
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end