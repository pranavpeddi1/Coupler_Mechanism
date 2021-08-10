%% Forward dynamic analysis of  Double crank
function fda()
global Minv h
load sckadata
clc;
m1 = 0.1;m2 = 0.2;m3 = 0.1;g = 9.81;
l1 = 0.4;l2 = 0.5;l3 =0.4;
J1 = m1*l1^2/12;J2 = m2*l2^2/12;J3 = m3*l3^2/12;
M = diag([m1 m1 J1 m2 m2 J2 m3 m3 J3]);
Minv = inv(M);
h = [0 -m1*g 0 0 -m2*g 0 0 -m3*g 0]';
Y0 = [pcoordsall(:,1)' zeros(1,9)]';
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
tspan = 0:0.001:20
[T,Y] = ode45(@eqm,tspan,Y0,options);
figure
plot(T,Y(:,3));% position of crank
title('position analysis of crank');
xlabel('time');
ylabel('position');
figure
plot(T,Y(:,9));% position of lever
title('position analysis of lever');
xlabel('time');
ylabel('position');
end

function dy = eqm(t,y)
global Minv h
[Phi,D] = constraints(t,y(1:9,1));
alpha = 0;beta = 0;
lambda = (D*Minv*D')\(gamma(y(1:9,1),y(10:18,1))-2*alpha*D*y(10:18,1)-beta^2*Phi-D*Minv*h);
dy(1:9,1) = y(10:18,1);
dy(10:18,1) = Minv*(h+D'*lambda);
norm(Phi)
end

% function [Phi,D] = constraints(t,q)
% l1 = 0.4;l2 = 0.5;l3 = 0.4;
% x1 = q(1);y1 = q(2);phi1 = q(3);
% x2 = q(4);y2 = q(5);phi2 = q(6);
% x3 = q(7);y3 = q(8);phi3 = q(9);
% 
% r1 = [x1 y1]';r2 = [x2 y2]';r3 = [x3 y3]';
% s_1_O = [-l1/2 0]';s_1_A = [l1/2 0]';
% s_2_A = [-l2/2 0]';s_2_B = [l2/2 0]';
% s_3_B = [-l3/2 0]';s_3_C = [l3/2 0]';   
% 
% S_1_O = A(phi1)*s_1_O;S_1_A = A(phi1)*s_1_A;
% S_2_A = A(phi2)*s_2_A;S_2_B = A(phi2)*s_2_B;
% S_3_B = A(phi3)*s_3_B;S_3_C=A(phi3)*s_3_C;
% 
% Phi = [r1+S_1_O;
%        r2+S_2_A-r1-S_1_A;
%        r3+S_3_B-r2-S_2_B;
%       [0.05 0]'-r3-S_3_C  
%        ];
% S_1_O_r = [0 -1;1 0]*S_1_O;S_1_A_r = [0 -1;1 0]*S_1_A; % Finding Rotated 
% S_2_A_r = [0 -1;1 0]*S_2_A;S_2_B_r = [0 -1;1 0]*S_2_B;
% S_3_B_r = [0 -1;1 0]*S_3_B;S_3_C_r = [0 -1;1 0]*S_3_C;
% 
% D = [eye(2) S_1_O_r zeros(2,6);
%     -eye(2) -S_1_A_r eye(2) S_2_A_r zeros(2,3);
%     zeros(2,3) -eye(2) -S_2_B_r eye(2) S_3_B_r ;
%     zeros(2,6) -eye(2) -S_3_C_r ;
%     ];
% end
% function output = gamma(q,qdot)
% l1 = 0.15;l2 = 0.16;l3 = 0.17;
% phi1 = q(3);phi2 = q(6);phi3 = q(9);
% phi1dot = qdot(3);phi2dot = qdot(6);phi3dot = qdot(9);
% s_1_O = [-l1/2 0]';s_1_A = [l1/2 0]';
% s_2_A = [-l2/2 0]';s_2_B = [l2/2 0]';
% s_3_B = [-l3/2 0]';s_3_C = [l3/2 0]';   
% 
% S_1_O = A(phi1)*s_1_O;S_1_A = A(phi1)*s_1_A;
% S_2_A = A(phi2)*s_2_A;S_2_B = A(phi2)*s_2_B;
% S_3_B = A(phi3)*s_3_B;S_3_C = A(phi3)*s_3_C;
% 
% output = [S_1_O*phi1dot^2;
%          -S_1_A*phi1dot^2+S_2_A*phi2dot^2;
%          -S_2_B*phi2dot^2+S_3_B*phi3dot^2;
%          -S_3_C*phi3dot^2      
%          ];
% end


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
      [0.2 0]'-r3-S_3_C  
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

function output = gamma(q,qdot)
l1 = 0.4;
l2 = 0.5;
l3 = 0.4;
phi1 = q(3);phi2 = q(6);phi3 = q(9);
phi1dot = qdot(3);phi2dot = qdot(6);phi3dot = qdot(9);
s_1_O = [-l1/2 0]';s_1_A = [l1/2 0]';
s_2_A = [-l2/2 0]';s_2_B = [l2/2 0]';
s_3_B = [-l3/2 0]';s_3_C = [l3/2 0]';   

S_1_O = A(phi1)*s_1_O;S_1_A = A(phi1)*s_1_A;
S_2_A = A(phi2)*s_2_A;S_2_B = A(phi2)*s_2_B;
S_3_B = A(phi3)*s_3_B;S_3_C = A(phi3)*s_3_C;

output = [S_1_O*phi1dot^2;   % gamma vector
         -S_1_A*phi1dot^2+S_2_A*phi2dot^2;
         -S_2_B*phi2dot^2+S_3_B*phi3dot^2;
         -S_3_C*phi3dot^2      
         ];
end

function output = A(phi)
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end