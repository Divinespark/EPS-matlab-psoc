
clc  
close all
clear all
Convergence_Tolerance=1e-6;  %Newton Raphson Convergence Tolerance
 disp(['Convergence Tolerance, Ep:',num2str(Convergence_Tolerance)]);  
 
%-----------------Node, tributary data---------------
% Node data
point = ...
    [...
        1	3	0        0.065      0     0     1.06      0       0
        2	1	1.82   	 0.64       0     0       1       0       0
      ];
% Matrix formed by node data: column 1 node number, column 2 node type:1 - PQ 2 - PV 3 - balance
% column 3 load active, column 4 load reactive, column 5 generation active, column 6 generation reactive
% column 7 voltage real, column 8 voltage imaginary
% column 9 node-to-ground conductance

% Branch data
zhilu = ...
   [...
   %bus_i	%bus_j	 % R	 % X	    % B       % k
   1        2        0     0.034 	  0       	1;
   ];
						
% bus_i denotes the first node of a branch.						
% bus_j denotes the end node of the branch.						
% R and X denote the branch resistance and reactance, respectively, and B denotes the dana.						
% k denotes the transformer ratio of the first node to the last node.	

x=length(point(:,1)); % Number of nodes
y=length(zhilu(:,1)); % Number of branches

%Matrix assignment of initial values.
TYPE=point(:,2);%Assign column 2 of the point matrix to TYPE, similar to the following
U=point(:,7);  %U is the nodal voltage matrix
a=point(:,8);  %a is the nodal voltage phase angle matrix ****** here should be in radians
P=point(:,5)-point(:,3);  %P is the nodal active power
Q=point(:,6)-point(:,4);  %Q is the nodal reactive power
I=zhilu(:,1); %I is the starting node numbering matrix
J=zhilu(:,2);  %J is the terminating node numbering matrix
Rij=zhilu(:,3); %R is the line resistance
Xij=zhilu(:,4); %X is the line reactance
Zij=Rij+j*Xij;  %Zij is the line impedance
Bij=zhilu(:,5);  %B0 is the n*1st order line to ground electroneutral value
KT=1./zhilu(:,6);  %K is the transformer ratio of order y*1 of the ij branch, if k=1 means no transformer, K=1 is the standard ratio, k not equal to 1 is the non-standard ratio
tic  %Program runtime starts
%Find the nodal derivative matrix Y
Y=zeros(x);   %Calculating the nodal conductance matrix
for m=1:y
    Y(I(m),J(m))=-1/( Zij(m)*(KT(m)) );
    Y(J(m),I(m))=Y(I(m),J(m));
    Y(I(m),I(m))=Y(I(m),I(m))+1/( Zij(m)*(KT(m))^2 )+j*Bij(m)/2;
    Y(J(m),J(m))=Y(J(m),J(m))+1/( Zij(m) )+j*Bij(m)/2;
end
for m=1:x  %--Calculation of conductivity to earth--
    Y(m,m)=Y(m,m)+1i*point(m,9);
end
G=real(Y);
%Find the B' matrix and its inverse matrix B1
B=imag(Y);%Find the imaginary part of Y, the nodal electro-nano matrix
 ph=find(TYPE(:,1)==3);%Find the balance node number
 BB=B;%BB matrix as an intermediate variable
 BB(:,ph)=[];%The corresponding column of the balance node number is empty and eliminated
 BB(ph,:)=[];%Balance node number corresponding row empty
 B1=BB;
 B1=inv(B1);%B1 matrix is obtained by inverting B1
 %Find B'' and its inverse matrix B2
phpv=find(TYPE(:,1)>1);%Find the number of the non-PQ node
 BB=B;     %The following is the same as finding B'
 BB(:,phpv)=[];
 BB(phpv,:)=[];
 B2=BB;    
 B2=inv(B2);%Find the B2 matrix
 
 %Calculate the active power imbalance at each node deltaPi
 k=0;   %k is the number of iterations
 kp=0;  %Calculate the convergence sign of the P-inequality measure deltaPi (0 means no convergence, 1 means convergence)
 kq=0;  %Calculate the convergence sign of the U-inequality measure deltaQi (0 means no convergence, 1 means convergence)
 notph=find(TYPE(:,1)<3);%Find the unbalanced node number
 deltaPi=zeros(x-1,1);%deltaPi is the (x-1)*1st order matrix,x is the number of nodes
 pq=find(TYPE(:,1)==1);%Find the PQ node number
 pqnum=size(B2);
 pqnum=pqnum(1);%Find the number of PQ nodes (since the dimension of the B1 matrix is equal to the number of PQ nodes)
 deltaQi=zeros(pqnum,1);%deltaQi is a pqnum*1st order matrix
 Convergence=[];
while((~kq|~kp)&(k<100)) %Cycle until 100 times
    k=k+1;
     for m=1:(x-1)%ÇódeltaPi
         sum1=0;
         for n=1:x
             sum1=sum1+U(notph(m))*U(n)*(G(notph(m),n)*cos(a(notph(m))-a(n))+B(notph(m),n)*sin(a(notph(m))-a(n)));
         end
         deltaPi(m)=P(notph(m))-sum1;
     end
     Up=U;   %Up is an intermediate variable
     Up(ph)=[];%Leave the row where the balance node is located empty
     Unotph=Up;%Find the vector of voltage columns except for the balance node
     deltaa=(-B1*(deltaPi./Unotph))./Unotph;%Find the inequality of the phase angle a
      for m=1:(x-1)   %Find the new iteration value matrix for the phase angle a
         a(notph(m))=a(notph(m))+deltaa(m);
     end
     max1=abs(deltaPi(1)/U(notph(1)));%Find the maximum value of the absolute value of deltaP/U
     for m=1:(x-2)
        % if abs(deltaPi(m)/U(notph(m)))<abs(deltaPi(m+1)/U(notph(m+1)))
          if max1<abs(deltaPi(m+1)/U(notph(m+1)))
             max1=abs(deltaPi(m+1)/U(notph(m+1)));
         end
     end
     if max1<=Convergence_Tolerance   %If the maximum value is satisfied, then kp is set to "1", indicating convergence
         kp=1;
     end
     for m=1:pqnum  %ÇódeltaQi
         sum2=0;
         for n=1:x
             sum2=sum2+U(pq(m))*U(n)*(G(pq(m),n)*sin(a(pq(m))-a(n))-B(pq(m),n)*cos(a(pq(m))-a(n)));
         end
         deltaQi(m)=Q(pq(m))-sum2;
     end
      Uq=U;%Uq is an intermediate variable
      Uq(phpv)=[];%Empty the row where the non-PQ node is located
      Upq=Uq;%Find the voltage column vector including the voltage at the PQ node
     deltaU=-B2*(deltaQi./Upq);%Find the inequality measure deltaU of U
     max2=max(abs(deltaQi./Upq)); %Find the maximum value of the absolute value of deltaQ/U
     if max2<=Convergence_Tolerance     %If the maximum value is satisfied, then kq is set to "1", indicating convergence
         kq=1;
     end
     for m=1:pqnum   %Find the new value of the iteration of U
         U(pq(m))=U(pq(m))+deltaU(m);
     end
     
     Convergence=[Convergence;max(max1,max2)];
     
 disp(['Number of iteration:',num2str(k),'  Max(|DP|,|DQ|):',num2str(max(max1,max2)),'  Voltage Magnitude at bus 2:',num2str(U(2,1)),'  Voltage Angle at bus 2:',num2str(a(2,1)*180/pi)]);
   
     
end

 %Find the line power Sij and Sji
 Sij=zeros(y,1);
 Sji=zeros(y,1);
 for m=1:y
     yi0=1/( Zij(m)*(KT(m))^2 )*(1-(KT(m)))+ i*Bij(m)/2;
     yj0=1/( Zij(m)*(KT(m)) )*((KT(m))-1)+ i*Bij(m)/2;
     Sij(m,1)=U(I(m))^2*( conj(yi0) ) + ( U(I(m))^2-U(I(m))*U(J(m))*( cos(a(I(m))-a(J(m)))+1i*sin(a(I(m))-a(J(m))) ) )*conj(-Y(I(m),J(m)));
     Sji(m,1)=U(J(m))^2*( conj(yj0) ) + ( U(J(m))^2-U(J(m))*U(I(m))*( cos(a(J(m))-a(I(m)))+1i*sin(a(J(m))-a(I(m))) ) )*conj(-Y(J(m),I(m)));
 end
 
deltaSij=Sij+Sji; %Find the power loss deltaSij
S=zeros(x,1)+i*zeros(x,1);%S is the nodal complex power, and is set to 0;
for b=1:x 
   for m=1:y
        if I(m)==b  %In the case of the starting node of a branch, this is calculated by the following formula
            S(b)=S(b)+Sij(m);
        elseif J(m)==b   %In the case of a branch termination node, this is calculated by the following equation
            S(b)=S(b)+Sji(m);   
        else  
            S(b)=S(b);
        end
  end
end
P=real(S);      %Nodal active power
Q=imag(S);      %Nodal reactive power

disp(['Active Power Injection at bus 2:',num2str(P(2,1)),'  Reactive Power Injection at bus 2:',num2str(Q(2,1))]);

figure(1);
kk=1:k;
plot(kk,Convergence(kk,1),'-b'),xlabel('k'),ylabel('Convergence Tolerance'); hold on;
