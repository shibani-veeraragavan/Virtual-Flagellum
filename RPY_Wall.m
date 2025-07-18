function Filament = RPY_Wall(Filaments,mu,wall)
% RPY  Solves the Stokes flow problem using the Rotne-Prager-Yamakawa
%      tensor.
%
%   RPY(Filaments,mu)
%   sets the velocities and angular velocities of the Filament objects
%   for given forces and torques, for spherical particles, i, of radius 
%   Filament.R(i) at positions Filament.X(i) in an unbounded Newtonian 
%   fluid of viscosity mu.
%
%   Details of RPY solver: Wajnryb et al., 2013 Journal of Fluid Mechanics,
%   "Generalization of the Rotne-Prager-Yamakawa mobility and shear
%   disturbance tensors".
%
%   These expressions correpond to equations (23)-(26) in the paper.


N = length(Filaments);

for i=1:N   % For each filament i
    
    Xi = Filaments(i).X;
    
    Ri = Filaments(i).R;
    
    V = zeros(3,Filaments(i).N_w);
    
    Omega = zeros(3,Filaments(i).N_w);
    
    for j=1:N   % For each filament j
        
        Xj = Filaments(j).X;
        
        Fj = Filaments(j).F(1:3,:);
        
        Tj = Filaments(j).F(4:6,:);
        
        Rj = Filaments(j).R;
        
        for m=1:Filaments(i).N_w
            
            Xm = Xi(:,m);
            
            am = Ri(m);
                     
            for n=1:Filaments(j).N_w
                
                Xn = Xj(:,n);
                
                R = ((Xm(1) - Xn(1))*(Xm(1) - Xn(1)) + (Xm(2) - Xn(2))*(Xm(2) - Xn(2)) + (Xm(3) - Xn(3))*(Xm(3) - Xn(3)))^0.5 + 10^-16;
                
                Xhat = (Xm - Xn); % At Xm due to Xn
                
                h = Xn(3); % Height above wall at z=0, positive half-space
                
                Xn_i = Xn - [0; 0; 2*h];    % Image point
                
                r = ((Xm(1) - Xn_i(1))*(Xm(1) - Xn_i(1)) + (Xm(2) - Xn_i(2))*(Xm(2) - Xn_i(2)) + (Xm(3) - Xn_i(3))*(Xm(3) - Xn_i(3)))^0.5 + 10^-16;
                
                xhat = (Xm - Xn_i);
                
                F = Fj(:,n);
                
                T = Tj(:,n);
                
                an = Rj(n);
               
                if m==n && i==j     % self-interaction
                    
                    V(:,m) = V(:,m) + F/(6*pi*mu*an);
                    Omega(:,m) = Omega(:,m) + T/(8*pi*mu*an^3);                    
                    if wall > 0
                        V(:,m) = V(:,m) + ((-9/16*an/h + 1/8*(an/h)^3 - 1/16*(an/h)^5)*[1 0 0; 0 1 0; 0 0 2] + 1/4*(an/h)^3*[0 0 0; 0 0 0; 0 0 1])*F/(6*pi*mu*an);
                        V(:,m) = V(:,m) - 1/(64*pi*mu)*(an^2/h^4)*[T(2); -T(1); 0];
                        Omega(:,m) = Omega(:,m) + 3*an^2/(32*h^4)*[F(2); -F(1); 0]/(6*pi*mu*an);
                        Omega(:,m) = Omega(:,m) + an^3/(16*h^3)*[-5 0 0; 0 -5 0; 0 0 -2]*T/(8*pi*mu*an^3);
                    end
                    
                else    % pair-wise interaction
                    
                    V(:,m) = V(:,m) + 1/(6*pi*mu*an)*((3/4*an/R + 1/2*an^3/R^3)*F + (3/4*an/R^3 - 3/2*an^3/R^5)*(F(1)*Xhat(1) + F(2)*Xhat(2) + F(3)*Xhat(3))*Xhat);
                    V(:,m) = V(:,m) - 1/(8*pi*mu*R^3)*(Xhat(1)*[0 0 0; 0 0 1; 0 -1 0] + Xhat(2)*[0 0 -1; 0 0 0; 1 0 0] + Xhat(3)*[0 1 0; -1 0 0; 0 0 0])*T;
                    Omega(:,m) = Omega(:,m) + 3/4*an/R^3*(Xhat(1)*[0 0 0; 0 0 1; 0 -1 0] + Xhat(2)*[0 0 -1; 0 0 0; 1 0 0] + Xhat(3)*[0 1 0; -1 0 0; 0 0 0])*F/(6*pi*mu*an);
                    Omega(:,m) = Omega(:,m) + an^3/2*(-1/R^3*[1 0 0; 0 1 0; 0 0 1] + 3/R^5*(Xhat*Xhat'))*T/(8*pi*mu*an^3);
                    
                    if wall > 0
                        MOF = (3/4*an/r^3*(xhat(1)*[0 0 0; 0 0 1; 0 -1 0] + xhat(2)*[0 0 -1; 0 0 0; 1 0 0] + xhat(3)*[0 1 0; -1 0 0; 0 0 0]) + (-3*an*h/(2*r^3) + 3*an^3*xhat(3)/(2*r^5))*[0 1 0; -1 0 0; 0 0 0]  ...
                            + (-9*an*h/(2*r^5) + 15*an^3*xhat(3)/(2*r^7))*([-xhat(2); xhat(1); 0]*xhat') + (9*an*h*xhat(3)/r^5 - 15*an^3*xhat(3)^2/r^7 + 3*an^3/(2*r^5))*([-xhat(2); xhat(1); 0]*[0,0,1]));
                        V(:,m) = V(:,m) + ((3/4*an*(-1/r - 2*h*xhat(3)/r^3 + 2*h^2/r^3) - an^3/2*(1/r^3 - 3*xhat(3)^2/r^5)- an^5/4*(-2/r^5 + 10*xhat(3)^2/r^7))*[1 0 0; 0 1 0; 0 0 1] ...
                            + (3/4*an*(-1/r^3 + 6*h*xhat(3)/r^5 - 6*h^2/r^5) - an^3/2*(-3/r^5 + 15*xhat(3)^2/r^7)- an^5/4*(-70*xhat(3)^2/r^9 + 10/r^7))*(xhat*xhat') ...
                            + (3/4*an*(2*h/r^3) - an^5/4*(20*xhat(3)/r^7))*([0; 0; 1]*xhat') ...
                            + (3/4*an*(2*h/r^3 - 12*h*xhat(3)^2/r^5 + 12*h^2*xhat(3)/r^5) - an^3/2*(6*xhat(3)/r^5 - 30*xhat(3)^3/r^7)- an^5/4*(140*xhat(3)^3/r^9 - 40*xhat(3)/r^7))*(xhat*[0, 0, 1]) ...
                            + (3/4*an*(-4*h^2/r^3) - an^3/2*(6*xhat(3)^2/r^5)- an^5/4*(8/r^5 - 60*xhat(3)^2/r^7))*[0 0 0; 0 0 0; 0 0 1])*F/(6*pi*mu*an);
                        V(:,m) = V(:,m) + MOF'*T/(8*pi*mu*an^3);
                        Omega(:,m) = Omega(:,m) + MOF*F/(6*pi*mu*an);
                        Omega(:,m) = Omega(:,m) + ((2/r^3 - 12*xhat(3)^2/r^5)*[1 0 0; 0 1 0; 0 0 1] - 6/r^5*(xhat*xhat') + 12/r^5*([-xhat(2); xhat(1); 0]*[-xhat(2); xhat(1); 0]') + 12*xhat(3)/r^5*([0;0;1]*xhat'))*T/(32*pi*mu);
                    end
                    
                end
                
            end
            
        end
        
    end
    
    Filaments(i).V(:,:) = [V;Omega];
    
end
Filament = Filaments(1);
end % End function.