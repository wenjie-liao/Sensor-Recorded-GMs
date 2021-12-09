function [ Ts,MaxDisp,MaxVel,MaxAcc ] = Spectrum_0519( A_acc )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

	% ***********linear response spectrum***********
	% initialization 
    Damp=0.05;            % damping ratio 0.05
	Dt=A_acc(2,1)-A_acc(1,1); % interval of acceleration records
	T_max = 6; % max period, unit:s
	dT_spect = 0.05; % interval of period, unit:s
	Ts = 0.0:dT_spect:T_max; % periods
    count = length(A_acc(:,1)); % acceleration length
    
    % Maximum relative displacement, velocity, maximum absolute acceleration
    MaxDisp = zeros(1,length(Ts));
    MaxVel = zeros(1,length(Ts));
    MaxAcc = zeros(1,length(Ts));
	
    t = 1;
	for T=0.0:dT_spect:T_max
		Frcy=2*pi/T ; % structural frequency
		DamFrcy=Frcy*sqrt(1-Damp*Damp);
		e_t=exp(-Damp*Frcy*Dt);
		s=sin(DamFrcy*Dt);
		c=cos(DamFrcy*Dt);
		A=zeros(2,2);             
		A(1,1)=e_t*(s*Damp/sqrt(1-Damp*Damp)+c);
		A(1,2)=e_t*s/DamFrcy;
		A(2,1)=-Frcy*e_t*s/sqrt(1-Damp*Damp);
		A(2,2)=e_t*(-s*Damp/sqrt(1-Damp*Damp)+c);
		d_f=(2*Damp^2-1)/(Frcy^2*Dt);
		d_3t=Damp/(Frcy^3*Dt);
		B=zeros(2,2);
		B(1,1)=e_t*((d_f+Damp/Frcy)*s/DamFrcy+(2*d_3t+1/Frcy^2)*c)-2*d_3t;
		B(1,2)=-e_t*(d_f*s/DamFrcy+2*d_3t*c)-1/Frcy^2+2*d_3t;
		B(2,1)=-e_t*(((Damp/(Frcy*Dt)+1)*s/DamFrcy)+(1/(Frcy^2*Dt))*c)+1/(Frcy^2*Dt);
		B(2,2)=e_t*((Damp/(Frcy*Dt)*s/DamFrcy)+(1/(Frcy^2*Dt))*c)-1/(Frcy^2*Dt); 
        
        Displace=zeros(1,count);  % relative displacement
        Velocity=zeros(1,count);  % relative velocity
        AbsAcce=zeros(1,count); % absolute acceleration
		for i=1:(count-1)
            Displace(i+1)=A(1,1)*Displace(i)+A(1,2)*Velocity(i)+B(1,1)*A_acc(i,2)+B(1,2)*A_acc(i,2);
			Velocity(i+1)=A(2,1)*Displace(i)+A(2,2)*Velocity(i)+B(2,1)*A_acc(i,2)+B(2,2)*A_acc(i,2);
			AbsAcce(i+1)=-2*Damp*Frcy*Velocity(i+1)-Frcy^2*Displace(i+1);
        end
        
        MaxDisp(1,t) = max(abs(Displace));
        MaxVel(1,t) = max(abs(Velocity));
        if T==0.0
            MaxAcc(1,t) = max(abs(A_acc(:,2)));
        else
            MaxAcc(1,t) = max(abs(AbsAcce));
        end
        t = t+1;
        
    end
    
end

