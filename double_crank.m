function double_crank
%% ------------------------------------------------------------------
%% Kinematic analysis of double crank mechanism in body-coordinates
%% Position, velocity, and acceleration analysis from t = 0 to t = 20 s
%% Coordinate-partitioning method
%% -------------------------------------------------------------------
clc
clear all
%% Position Analysis
q = [0 0.2 pi/2 0.25 0.5 0 0.2 0.25 pi/2]';% variables x1,y1,phi1,x2,y2,phi2,...phi3.
t = 0:0.001:0.01;% time
N = size(t,2);% Assigning length of time in 1*1 matrix form
pcoordsall = zeros(9,N);% storing variabales at each time step i.e, from 0 to 20s in the intervals of 0.001

for i = 1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i)); %driver constraints ,t(i)=time
    q(3) = phi1;   
   [pcoordsall(:,i),q] = NewtonRaphson(q,9,1e-10,120,@constraints,t(i)); % storing all variables for each time interval in pcoords. All variables are output got from newton raphson
   [Phi,D] = constraints(t(i),pcoordsall(:,i));
   norm(Phi)
end

figure;
plot(t,pcoordsall(3,:),'k')
title('position analysis');
xlabel('time');
ylabel('position');
hold on
plot(t,pcoordsall(9,:),'r')
legend('phi1','phi3');

%% Velocity Analysis
qdot = zeros(9,1); %velocity vector
vcoordsall = zeros(9,N);%storing velocities at each time step
for i = 1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i));
    qdot(3) = phi1dot;
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    Dnew = [D(:,1:2) D(:,4:9)];% third column is removed because we know the values of third column after multiplying with qdot(3)
    Rhs = -D(:,3)*qdot(3);% The removed third  column is moved to RHS of the Equation .
    qdotnew = Dnew\Rhs; % Unknown velocities are found by this . qdotnew has (8*1)
    vcoordsall(:,i) = [qdotnew(1:2,1)' qdot(3) qdotnew(3:8,1)']'; % Arranging qdot(3) to make velocity vector complte size of 9*1.
end
figure;
plot(t,vcoordsall(3,:),'r')
title('velocity analysis');
xlabel('time');
ylabel('velocity');
hold on
plot(t,vcoordsall(9,:),'k')
legend('phi1dot','phi3dot');

%%  Acceleration Analysis
qddot = zeros(9,1);
acoordsall = zeros(9,N);
for i = 1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i));
    qddot(3) = phi1ddot;
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    Dnew = [D(:,1:2) D(:,4:9)];% third column is removed because we know the values of third column after multiplying with qddot(3)
    rhs = -D(:,3)*qddot(3)+gamma(pcoordsall(:,i),vcoordsall(:,i)); % The removed third column is moved to RHS of the Equation .
    qddotnew = Dnew\rhs;
    acoordsall(:,i) = [qddotnew(1:2,1)' qddot(3) qddotnew(3:8,1)']';
end
figure;
plot(t,acoordsall(3,:),'b')
title('Acceleration analysis');
xlabel('time');
ylabel('Acceleratioin');
hold on
plot(t,acoordsall(9,:),'k')
legend('phi1ddot','phi3ddot');
save sckadata.mat t pcoordsall vcoordsall acoordsall %saving all variables data
end

%% right hand term of acceleration = Gamma
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

% Driver function

function [phi1,phi1dot,phi1ddot] = driver(t)
if t<=1
    phi1 = pi/2+t^3-t^4/2;
    phi1dot = 3*t^2-2*t^3;
    phi1ddot = 6*t-6*t^2;
else
    phi1 = pi/2+1-0.5+1*(t-1);
    phi1dot = 1;
    phi1ddot = 0;
end
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

function output = A(phi) % Defining Rotation matrix A
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end