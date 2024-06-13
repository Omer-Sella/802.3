function [return_struct]= capture_RIL_RILN(chdata)
% History:
% 1. 12th April, 2019 (Intial release)
%
% 2. 11th October, 2021
%   - Details:
%       1] Revised the equation of RIL for missing conj(rho_port1) and conj(rho_port2).
%       2] Revised the selection criteria for the solution of the quadratic
%       equation in finding the reflection coefficient (rho).
%   - Impact:
%       => Zero impact in |RIL|, while impact on angle(RIL).
%   - Previous:
%       %---start. For passive networks the reflection coefficient should be less than one
%       if all(abs(solution_1(idx_start:end))< 1) && ~all(abs(solution_2(idx_start:end))< 1)
%           rho_port2= solution_1;
%        elseif all(abs(solution_2(idx_start:end))< 1) && ~all(abs(solution_1(idx_start:end))< 1)
%          rho_port2= solution_2;
%        else
%            rho_port2= solution_1;
%           %     error('Please contact the tool developer. It appears a odd case has appeared.');
%       end
%       %---end. For passive networks the reflection coefficient should be less than one
%
%       RIL= conj((1- rho_port1).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- rho_port2).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );
%   - Change:
%       %---start. Given the real part of the impedance is to be positive
%       Z_solution_1= (solution_1*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_1);
%       Z_solution_2= (solution_2*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_2);
% 
%       rho_port2= zeros(length(solution_1), 1);
%        for solution_idx= 1:length(solution_1)
%          if real(Z_solution_1(solution_idx))>0 && real(Z_solution_2(solution_idx))<=0
%              rho_port2(solution_idx, 1)= solution_1(solution_idx);
%           elseif real(Z_solution_2(solution_idx))>0 && real(Z_solution_1(solution_idx))<=0
%              rho_port2(solution_idx, 1)= solution_2(solution_idx);
%          else
%               error('An odd case has occured. Please contact the tool developer.');        
%           end
%       end
%       %---end. Given the real part of the impedance is to be positive
%       RIL= conj((1- conj(rho_port1)).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- conj(rho_port2)).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );

% Definition:
% This function captures the reflectionless insertion loss (RIL) and reflective insertion loss nois (RILN) for any arbitary S-parameter

% Author:
% Hansel Dsilva (dsilvahansel@gmail.com or hanseldsilva@achronix.com)
% Acknowledgement: Adam Gregory, Richard Mellitz, Beomtaek Lee and Amit Kumar.

% This function has been shared by Hansel for others to evaluate for reflectionless insertion loss (RIL) and reflective insertion loss nois (RILN) for any arbitary S-parameter.

% Reference:
% 1] H. Dsilva et al., "Finding Reflective Insertion Loss Noise and Reflectionless Insertion Loss," 2020 DesignCon.
% 2] H. Dsilva, A. Jain, J. Sasikala and A. Kumar, "Novel Signal Integrity Application of Power Wave Scattering Matrix theory," 2019 IEEE MTT-S International Microwave and RF Conference (IMARC), 2019, pp. 1-7, doi: 10.1109/IMaRC45935.2019.9118701.

% Input:
% 1] SCH: S-matrix structure
% SCH.Frequencies= faxis;
% SCH.Parameters(1,1,:)= sdd11;
% SCH.Parameters(2,2,:)= sdd22;
% SCH.Parameters(1,2,:)= sdd12;
% SCH.Parameters(2,1,:)= sdd21;
% SCH.NumPorts= 2;
% SCH.Impedance= 100;

% Output: Struct returned has the following,
% return_struct.RIL          %Reflectionless Insertion Loss as a complex number
% return_struct.RIL_dB          %Reflectionless Insertion Loss in decibel
% return_struct.RILN          %Reflective Insertion Loss Noise as a complex number
% return_struct.RILN_dB          %Reflective Insertion Loss Noise in decibel
% return_struct.Z_port1          %Frequency dependent, complex values of impedance for termination of port 1
% return_struct.Z_port2          %Frequency dependent, complex values of impedance for termination of port 2
% return_struct.rho_port1          %Frequency dependent, complex values of the reflection coefficient of port 1
% return_struct.rho_port2          %Frequency dependent, complex values of reflection coefficient of port 2
% return_struct.freq          %Frequency axis

SCH.Parameters(1,1,:)=chdata(1).sdd11_orig;
SCH.Parameters(2,2,:)=chdata(1).sdd22_orig;
SCH.Parameters(1,2,:)=chdata(1).sdd12_orig;
SCH.Parameters(2,1,:)=chdata(1).sdd21_orig;
SCH.Frequencies=chdata(1).faxis;
SCH.Impedance=100;%<------------------in a future release, may want to parameterize this based on the read .sNp
SCH.NumPorts= 2;

%---start. allowed is only a 2-port network having a transmitter and receiver
if size(SCH.Parameters, 1)~=2
	fprintf('Size of given S matrix is %d.', size(SCH.Parameters, 3) );
	error('Allowed is only a 2-port network having a transmitter and receiver.');
end
%---end. allowed is only a 2-port network having a transmitter and receiver

%---start. do not include the DC point given sinusoidals at DC are not
%defined
if SCH.Frequencies(1)==0
    idx_start= 2;
else
    idx_start= 1;
end
%---end. do not include the DC point given sinusoidals at DC are not
%defined

Sdd11= squeeze(SCH.Parameters(1,1,idx_start:end));
Sdd12= squeeze(SCH.Parameters(1,2,idx_start:end));
Sdd21= squeeze(SCH.Parameters(2,1,idx_start:end));
Sdd22= squeeze(SCH.Parameters(2,2,idx_start:end));

a= -Sdd22+ Sdd11.*Sdd22.*conj(Sdd11)- Sdd21.*Sdd12.*conj(Sdd11);
b= 1+ Sdd22.*conj(Sdd22)+...
    Sdd11.*Sdd22.*conj(Sdd21).*conj(Sdd12)-...
    Sdd21.*Sdd12.*conj(Sdd21).*conj(Sdd12)-...
    Sdd11.*Sdd22.*conj(Sdd11).*conj(Sdd22)+...
    Sdd12.*Sdd21.*conj(Sdd11).*conj(Sdd22)-...
    Sdd11.*conj(Sdd11);
c= -conj(Sdd22)-...
    Sdd11.*conj(Sdd21).*conj(Sdd12)+...
    Sdd11.*conj(Sdd11).*conj(Sdd22);

solution_1= (-b+sqrt((b.^2)-(4*a.*c)))./(2*a);
solution_2= (-b-sqrt((b.^2)-(4*a.*c)))./(2*a);

clear a b c

%---start. Given the real part of the impedance is to be positive
Z_solution_1= (solution_1*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_1);
Z_solution_2= (solution_2*conj(SCH.Impedance)+ SCH.Impedance)./(1-solution_2);

rho_port2= zeros(length(solution_1), 1);
for solution_idx= 1:length(solution_1)
    if real(Z_solution_1(solution_idx))>0 && real(Z_solution_2(solution_idx))<=0
        rho_port2(solution_idx, 1)= solution_1(solution_idx);
    elseif real(Z_solution_2(solution_idx))>0 && real(Z_solution_1(solution_idx))<=0
        rho_port2(solution_idx, 1)= solution_2(solution_idx);
    else
        error('An odd case has occured. Please contact the tool developer.');        
    end
end
%---end. Given the real part of the impedance is to be positive

rho_port1= conj(Sdd11+ ((rho_port2.*Sdd21.*Sdd12)./(1-(rho_port2.*Sdd22))));

%---start. calculate the equivalent port impedance
Z_port1= (rho_port1*conj(SCH.Impedance)+ SCH.Impedance)./(1-rho_port1);
Z_port2= (rho_port2*conj(SCH.Impedance)+ SCH.Impedance)./(1-rho_port2);
%---end. calculate the equivalent port impedance


% %---start. The reflectionless insertion loss is the insertion loss corresponding
% %to zero reflections.
RIL= conj((1- conj(rho_port1)).*(sqrt(abs(1- rho_port1.*conj(rho_port1)))./abs(1-rho_port1)  )).*(  Sdd21.*(1- rho_port2.*Sdd22)+ (Sdd22- conj(rho_port2)).*rho_port2.*Sdd21  )./ ( ((1- conj(rho_port2)).*(sqrt(abs(1- rho_port2.*conj(rho_port2)))./abs(1-rho_port2)  )) .* (  -rho_port1.*rho_port2.*Sdd21.*Sdd12+  (1- rho_port2.*Sdd22).*(1- rho_port1.*Sdd11)  ) );
%---end. The reflectionless insertion loss is the insertion loss corresponding
%to zero reflections.

%---start. Calculate RILN (Reflective Insertion Loss Noise)
RILN= RIL- Sdd21;
RILN_dB= 20*log10(abs(Sdd21))- 20*log10(abs(RIL));
%---end. Calculate RILN (Reflective Insertion Loss Noise)

%---start.  preparing returns struct
return_struct.RIL= RIL;		%Reflectionless Insertion Loss as a complex number
return_struct.RIL_dB= 20*log10(abs(RIL));		%Reflectionless Insertion Loss in decibel
return_struct.RILN= RILN;		%Reflective Insertion Loss Noise as a complex number
return_struct.RILN_dB= RILN_dB;		%Reflective Insertion Loss Noise in decibel
return_struct.Z_port1= Z_port1;		%Frequency dependent, complex values of impedance for termination of port 1
return_struct.Z_port2= Z_port2;		%Frequency dependent, complex values of impedance for termination of port 2
return_struct.rho_port1= rho_port1;		%Frequency dependent, complex values of the reflection coefficient of port 1
return_struct.rho_port2= rho_port2;		%Frequency dependent, complex values of reflection coefficient of port 2
return_struct.freq= SCH.Frequencies(idx_start:end);		%Frequency axis
%---end. preparing returns struct