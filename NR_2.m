%---------Initial settings-----------
clc;
clear;
close all 
Convergence_Tolerance=1e-6;  %Newton Raphon Convergence Tolerance
 disp(['Convergence Tolerance, Ep:',num2str(Convergence_Tolerance)]);  

 %-----------------Node, tributary data---------------
% Node data
node = ...
    [...
        1	3	0        0.065      0     0     1.06      0       0
        2	1	1.82   	 0.64       0     0       1       0       0
      ];
% Matrix formed by node data: column 1 node number, column 2 node type:1 - PQ 2 - PV 3 - balance
% column 3 load active, column 4 load reactive, column 5 generation active, column 6 generation reactive
% column 7 voltage real, column 8 voltage imaginary
% column 9 node-to-ground conductance

% Branch data
branch = ...
   [...
 %bus_i   %bus_j    % R     % X	    % Bi/2    	% k
   1        2        0     0.034 	  0       	1;
   ];
						
% bus_i denotes the first node of a branch.						
% bus_j denotes the end node of the branch.						
% R and X denote the branch resistance and reactance, respectively, and B denotes the dana.						
% k denotes the transformer ratio of the first node to the last node.		

n=size(node,1);  %Number of nodes
l=size(branch,1);  %Number of branches

%%--------------------------------------
%Calculating the nodal derivative matrix
%%--------------------------------------
Y=zeros(n);   %Calculating the nodal conductance matrix
for j=1:l                      
    s1=branch(j,1);s2=branch(j,2);  %s1 is the first segment node of each branch and s2 is the end node
    kij=branch(j,6); %Branch transformer ratio
    Y(s1,s2)=-1/( (branch(j,3)+1i*branch(j,4))*kij );
    Y(s2,s1)=Y(s1,s2);
    Y(s1,s1)=Y(s1,s1)+1/( (branch(j,3)+1i*branch(j,4))*kij^2 )+i*branch(j,5);  
    Y(s2,s2)=Y(s2,s2)+1/(branch(j,3)+1i*branch(j,4))+i*branch(j,5);  
end
for a=1:n  %--Calculation of conductivity to earth--
    Y(a,a)=Y(a,a)+1i*node(a,9);
end
G=real(Y);   %Nodal derivative matrix real part matrix
B=imag(Y);   %Nodal derivative matrix imaginary part matrix


%%Newton Raphon method of calculating currents

e=node(:,7);    %Initial value of the real part of the voltage
f=zeros(n,1);   %Initial value of the imaginary part of the voltage

for a=1:n  %Search to balance nodes
    if node(a,2)==3
        slack=a;
    end
end

it=0;   %Number of iterations
Convergence=[];
while 1
    
   DP=zeros(n,1); DQ=zeros(n,1); %Uneven measurement
   Pi=zeros(n,1); Qi=zeros(n,1); %Nodal injected power calculated by the tidal equation
   for a=1:n   %Calculating Pi, Qi
       if node(a,2)==1 || node(a,2)==3 %Determined to be a PQ node
          for b=1:n
              Pi(a,1)=Pi(a,1)+e(a,1)*( G(a,b)*e(b,1)-B(a,b)*f(b,1) )+ f(a,1)*( G(a,b)*f(b,1)+B(a,b)*e(b,1) );
              Qi(a,1)=Qi(a,1)+f(a,1)*( G(a,b)*e(b,1)-B(a,b)*f(b,1) )- e(a,1)*( G(a,b)*f(b,1)+B(a,b)*e(b,1) );
          end
          DP(a,1)=( node(a,5)-node(a,3) )-Pi(a,1);  %Calculating unevenness
          DQ(a,1)=( node(a,6)-node(a,4) )-Qi(a,1);
       elseif node(a,2)==2  %Determined to be a PV node
           for b=1:n
              Pi(a,1)=Pi(a,1)+e(a,1)*( G(a,b)*e(b,1)-B(a,b)*f(b,1) )+ f(a,1)*( G(a,b)*f(b,1)+B(a,b)*e(b,1) );
              Qi(a,1)=e(a,1)^2+f(a,1)^2;
           end
          DP(a,1)=( node(a,5)-node(a,3) )-Pi(a,1);
          DQ(a,1)=node(a,7)^2-Qi(a,1);  %PV node voltage imbalance
       end  %Balance nodes regardless of
   end

   DPQ=zeros(2*n,1); %Integration of active, reactive and voltage inequalities into one matrix
   for a=1:n
       DPQ(2*a-1,1)=DP(a,1);
       DPQ(2*a,1)=DQ(a,1);
   end
   DPQ([2*slack-1,2*slack],:)=[];  % The slack node is a balancing node and does not participate in iterations, remove
   
    xV=zeros(n,1);   %Node voltage amplitude
    xR=zeros(n,1);   %Node voltage phase angle
    for a=1:n
       xV(a,1)=sqrt( e(a,1)^2+f(a,1)^2 );
       xR(a,1)=atan(f(a,1)/e(a,1))*180/pi;
    end
   disp(['Number of iteration:',num2str(it),'  Max(|DP|,|DQ|):',num2str(max(abs(DPQ))),'  Voltage Magnitude at bus 2:',num2str(xV(2,1)),'  Voltage Angle at bus 2:',num2str(xR(2,1))]);
   
   %If the inequality is less than the convergence accuracy, then the convergence condition has been reached and the iteration is skipped
   Convergence=[Convergence;max(abs(DPQ))];
   if max(abs(DPQ))<=Convergence_Tolerance
       break;
   end
   
   %If accuracy is achieved, the Jacobi matrix and corrections are calculated downwards and the iteration continues 
    H=zeros(n,n);  % Calculation of DP/De
    N=zeros(n,n);  % Calculation of DP/Df
    K=zeros(n,n);  % Calculate DQ/De for the PQ node and DU/De for the PV node
    L=zeros(n,n);  % Calculate DQ/Df for the PQ node and DU/Df for the PV node
    for a=1:n
       if node(a,2)==1 %Determined to be a PQ node
           for b=1:n
               if a~=b 
                   H(a,b)=-( G(a,b)*e(a,1)+B(a,b)*f(a,1) );
                   N(a,b)=B(a,b)*e(a,1)-G(a,b)*f(a,1) ;
                   K(a,b)=B(a,b)*e(a,1)-G(a,b)*f(a,1) ;
                   L(a,b)=G(a,b)*e(a,1)+B(a,b)*f(a,1) ;
               elseif a==b
                   for c=1:n
                      H(a,b)=H(a,b)-( G(a,c)*e(c,1)-B(a,c)*f(c,1) ); 
                      N(a,b)=N(a,b)-( G(a,c)*f(c,1)+B(a,c)*e(c,1) );
                      K(a,b)=K(a,b)+( G(a,c)*f(c,1)+B(a,c)*e(c,1) );
                      L(a,b)=L(a,b)-( G(a,c)*e(c,1)-B(a,c)*f(c,1) );
                   end
                   H(a,b)=H(a,b)-G(a,a)*e(a,1)-B(a,a)*f(a,1);
                   N(a,b)=N(a,b)+B(a,a)*e(a,1)-G(a,a)*f(a,1);
                   K(a,b)=K(a,b)+B(a,a)*e(a,1)-G(a,a)*f(a,1);
                   L(a,b)=L(a,b)+G(a,a)*e(a,1)+B(a,a)*f(a,1);
               end
           end
       elseif node(a,2)==2 %Determined to be a PV node
           for b=1:n
               if a~=b 
                   H(a,b)=-( G(a,b)*e(a,1)+B(a,b)*f(a,1) );
                   N(a,b)=B(a,b)*e(a,1)-G(a,b)*f(a,1) ;
                   K(a,b)=0 ; 
                   L(a,b)=0 ;
               elseif a==b
                   for c=1:n
                      H(a,b)=H(a,b)-( G(a,c)*e(c,1)-B(a,c)*f(c,1) ); 
                      N(a,b)=N(a,b)-( G(a,c)*f(c,1)+B(a,c)*e(c,1) );
                   end
                   H(a,b)=H(a,b)-G(a,a)*e(a,1)-B(a,a)*f(a,1);
                   N(a,b)=N(a,b)+B(a,a)*e(a,1)-G(a,a)*f(a,1);
                   K(a,b)=-2*e(a,1); 
                   L(a,b)=-2*f(a,1);
               end
           end
       end  %Balance nodes regardless of
    end
    jacobi=zeros(2*n); %Jacobi Matrix
    for a=1:n
        for b=1:n
           jacobi(2*a-1,2*b-1)=H(a,b);
           jacobi(2*a-1,2*b)=N(a,b);
           jacobi(2*a,2*b-1)=K(a,b);
           jacobi(2*a,2*b)=L(a,b);
        end
    end
    jacobi([2*slack-1,2*slack],:)=[]; jacobi(:,[2*slack-1,2*slack])=[]; % The slack node is a balancing node and does not participate in iterations, remove
    
    delt=-jacobi\DPQ;    %Calculation of corrections to variables
    
    %Balance node (node with slack number) with e and f corrections of 0
    if slack==1
        delt=[0;0;delt];
    elseif slack==n
        delt=[delt;0;0];
    else
        delt1=delt(1:2*slack-2,:);
        delt2=delt(2*slack-1:2*n-2,:);
        delt=[delt1;0;0;delt2];
    end
      
    
    for a=1:n
        e(a,1)=e(a,1)+delt(2*a-1,1);
        f(a,1)=f(a,1)+delt(2*a,1);
    end

      
    it=it+1; %Number of iterations
   
end

Pi=zeros(n,1); Qi=zeros(n,1);
for a=1:n   %Calculating Pi, Qi
    for b=1:n
        Pi(a,1)=Pi(a,1)+e(a,1)*( G(a,b)*e(b,1)-B(a,b)*f(b,1) )+ f(a,1)*( G(a,b)*f(b,1)+B(a,b)*e(b,1) );
        Qi(a,1)=Qi(a,1)+f(a,1)*( G(a,b)*e(b,1)-B(a,b)*f(b,1) )- e(a,1)*( G(a,b)*f(b,1)+B(a,b)*e(b,1) );
    end
end

disp(['Active Power Injection at bus 2:',num2str(Pi(2,1)),'  Reactive Power Injection at bus 2:',num2str(Qi(2,1))]);

figure(1);
k=1:it+1;
plot(k-1,Convergence(k,1),'-b'),xlabel('k'),ylabel('Convergence Tolerance'); hold on;


