% n=Body coordinates , For double crank it is 9, q= variables =9 ,t=time vector
function [sol,q] = NewtonRaphson(q,n,tol,iter_max,constraints,t)
coords = zeros(iter_max,n+1);% here n+1 is used because the columns includes error corresponding to the variables q with respective to ith iteration
flag = 0;sol = [];
for i = 1:iter_max
   [Phi,D] = constraints(t,q);   
   D = [D(:,1:2) D(:,4:9)];% 3 column is eliminated because phi1 is driver constraint so body 1 moves with given i/p and now D is 8*8
   err = sqrt(Phi'*Phi);
   coords(i,:) = [err q'];
   if err<tol
       flag = 1;
       sol = coords(i,2:n+1)';
       break;
   end
   delta_q = -D\Phi;% Back slash inverse of D
   delta_q = [delta_q(1:2,1)' 0 delta_q(3:8,1)']'; % Since phi1 is removed so to make 9*1 size of variable matrix 0 is added
   q = q + delta_q;
end
if flag == 0
    'Convergence failed in Newton-Raphson'
    return;
end
end