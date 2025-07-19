function Filament = Ani_Wall(Filament,mu,wall)
% Computes the translational and angular velocities of the filament's segments given the forces and torques 
% using Resistive Force Theory, with coefficients as provided by Lighthill, SIAM Review 18, 161–230 (1976) 
% and Katz et al., J Fluid Mech 72, 529 (1975). 

L = Filament.Length;
dL = Filament.DL;
Xi = Filament.X;
Ri = Filament.R;
qi = Filament.Q;

Fi = Filament.F(1:3,:);
Ti = Filament.F(4:6,:);

% Initialise velocity vectors
V = zeros(3,Filament.N_w);
Omega = zeros(3,Filament.N_w);
        
for m=1:Filament.N_w

    a = Ri(m);
    Xn = Xi(:,m);
    F = Fi(:,m);
    T = Ti(:,m);

    if wall>0
        % Use wall coefficients of Katz et al.
        h = Xn(3); % Height above wall at z=0, positive half-space
        if h/L>10
            Ct_inv = log(2*L/a)/(2*pi*mu*dL);
            Cn_inv = (log(2*L/a)+0.5)/(4*pi*mu*dL);
        elseif h/L<=10 && h/L>1
            Ct_inv = (log(2/a)-0.807-3*L/8/h)/(2*pi*mu*dL);
            Cn_inv = (log(2/a)+0.193-3*L/4/h)/(4*pi*mu*dL);
        else
            Ct_inv = log(2*h/a)/(2*pi*mu*dL);
            Cn_inv = log(2*h/a)/(4*pi*mu*dL);
        end
    else
        % Use free-space coefficients of Lighthill
        Ct_inv = log(0.09*L/a)/(2*pi*mu*dL);
        Cn_inv = (log(0.09*L/a)+0.5)/(4*pi*mu*dL);
    end

    qm = qi(m,1:4);
    d1m = QuaternionRotation(qm,[1;0;0]);
    d2m = QuaternionRotation(qm,[0;1;0]);
    d3m = QuaternionRotation(qm,[0;0;1]);
    
    Fvec = [sum(F.*d1m)*Ct_inv; sum(F.*d2m)*Cn_inv; sum(F.*d3m)*Cn_inv];
    V(1:3,m) = QuaternionRotation(qm, Fvec);

    Omega(1,m) = Omega(1,m) + T(1)/(8*pi*mu*amax^3);
    Omega(2,m) = Omega(2,m) + T(2)/(8*pi*mu*amax^3);
    Omega(3,m) = Omega(3,m) + T(3)/(8*pi*mu*amax^3);

end

Filament.V(:,:) = [V;Omega];

end
