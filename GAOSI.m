%====================================================================
% Gauss Seidel Power Flow Solution
% Program code for power flow solution for a two-bus power network
% Bus 1: Slack bus, V and Theta are known 
% Bus 2: PQ bus, V and Theta need to be determined
% Written by Prof. Xiao-Ping Zhang
%====================================================================

% Step 1: For Bus Admittance Matrix 
z12 = 0.0+0.034i;
y12 = 1/z12;
Y11 = y12;
Y12 = -y12;
Y21 = Y12;
Y22 = y12;

% Step 2: Power Injections: The number in the bracket indicates bus number
PQ = [];
% Load at bus 2 is 0.8 +j0.25; so the power injection is:
PQ(1) = -0.0-0.065i
PQ(2) = -1.82-0.64i

% Step 3: Initilization of power flow soultions
VK = [];
VK1 = [];
Ik = [];
% Voltage at bus 1 (slack bus) is given as:
VK(1) = 1.06+0.0i;
VK1(1) = 1.06+0.0i;
% Voltage at bus 2 (flat start)
VK(2) = 1.0+0.0i
VK1(2) = 1.0+0.0i
IK(2)= 0.0+0.0i

%Step 4: Enter convergence tolerance
Ep=input('Convergence Tolerance, Ep:');
fprintf('\n');
fprintf('%12.5e\n',Ep);

%Step 5: Gauss Seidel Iterations
for k=1:1000

% Step 5.1: IK(2) = conjugate ((P(2) + jQ(2))/VK(2)) 
IK(2) = conj(PQ(2)/VK(2)); 

% Step 5.2: Voltage at Bus 2 (PQ bus) VK+1(2) = IK(2) - Y21xV1
VK1(2) = (IK(2) - Y21*VK1(1))/Y22;
fprintf('Number of iteration: %4i Max |Vk+1 - Vk|= %12.5e Voltage Magnitude at bus 2: %12.9e Voltage Angle at bus 2: %12.9e\n', k, abs(VK1(2)-VK(2)), abs(VK1(2)), angle(VK1(2))*180/pi);

% Step 5.3: Determine the Convergence
  if (abs(VK1(2)-VK(2))< Ep)
 % If max |Vk+1 - Vk| < Ep, stop
  break;
 else
 % Otherwise continue the iteration and update solution: Vk = Vk+1
 VK(2) = VK1(2);
 end;

end

% Step 6: Calculate current injcetion at Slack bus = Bus 1:
IK(1)= Y11*VK1(1) + Y12*VK1(2);

% Step 7: Calculate power injcetion at Slack bus = Bus 1:
PQ(1) = VK1(1)*conj(IK(1)); 

% Step 8: Output power injections at Bus 1 (Slack bus)
fprintf('Active Power Injection at bus 2: %12.5e Reactive Power Injection at bus 2: %12.5e\n', real(PQ(1)), imag(PQ(1)));