function [E] = J_block_E(Filament,dt,mu)
% J_BLOCK_E(Filament,dt,mu)   Provides the terms relating the Lie algebra 
%                             update equations on the Lie algebra elements.
% Rows correspond to m=1...N-1, Columns correspond to n=1...N-1.
% Input: Filament (object), dt (scalar), mu (scalar).
% Output: C_mn = df^2_m/du_n (Size 3N x 3N) 

N_w = Filament.N_w;

E = zeros(3*N_w);

DL = Filament.DL;

KB = Filament.KB + 3*Filament.NB/dt;

KT = Filament.KT + 3*Filament.NT/dt;

strain_twist = Filament.StrainTwist;

Q = Filament.Q;

Lambda = Filament.Lambda;

U = Filament.U;

Omega = Filament.V(4:6,:);

R = Filament.R;

% ma = Filament.ma; %

% We deal with the end-particle cases seperately, because we need D1_m-1,
% D1_m and D1_m+1 of which only two are present at the ends.

% n = 1

Tfac = 1/(8*pi*mu*R(1)^3);

D1m = QuaternionRotation(Q(1,1:4),[1;0;0]);
D2m = QuaternionRotation(Q(1,1:4),[0;1;0]);
D3m = QuaternionRotation(Q(1,1:4),[0;0;1]);

D1m1 = QuaternionRotation(Q(2,1:4),[1;0;0]);
D2m1 = QuaternionRotation(Q(2,1:4),[0;1;0]);
D3m1 = QuaternionRotation(Q(2,1:4),[0;0;1]);

% Diagonal block

constraint_part = -0.5*DL*rcross(Lambda(:,1))*rcross(D1m);
% constraint_part = 0.5*DL*rcross(Lambda(:,1))*rcross(D1m); %

beta = -strain_twist(3,1) + 0.5*(D2m1' * D3m - D2m' * D3m1)/DL;

elastic_part = 0.5*(D1m + D1m1)*(cross(D3m,D2m1) - cross(D2m,D3m1))'/DL;
 
elastic_part = 0.5*KT*(beta*rcross(D1m) + elastic_part);

elastic_part = elastic_part + KB*(rcross(D1m1)*rcross(D1m)/DL...
    - 0.5*strain_twist(1,1)*rcross(D2m)...
    - 0.5*strain_twist(2,1)*rcross(D3m));

dTm_dun = Tfac*(elastic_part + constraint_part);

Lie_mat = rcross(U(1:3,1))';    % [u_n x] lie algebra
Omega_mat = rcross(Omega(:,1)); % [x w_n] angular velocity

block = -2*dt*(dTm_dun - 0.5*(Omega_mat + Lie_mat*dTm_dun) + ...
    (rcross(cross(U(1:3,1),Omega(:,1))) + Lie_mat*Omega_mat + Lie_mat^2*dTm_dun)/12)/3;

E(1:3,1:3) = block;

% Off-diagonal block

elastic_part = -0.5*(D1m + D1m1)*(cross(D3m,D2m1) - cross(D2m,D3m1))'/DL;

elastic_part = 0.5*KT*(beta*rcross(D1m1) + elastic_part);

elastic_part = elastic_part + KB*(rcross(D1m)'*rcross(D1m1)/DL...
    - 0.5*strain_twist(1,1)*rcross(D2m1)...
    - 0.5*strain_twist(2,1)*rcross(D3m1));

dTm_dun = Tfac*elastic_part;

block = -2*dt*(dTm_dun - 0.5*Lie_mat*dTm_dun + Lie_mat^2*dTm_dun/12)/3;

E(1:3,4:6) = block; % Placed at n=m+1

% n = N_w

Tfac = 1/(8*pi*mu*R(N_w)^3);

D1m0 = QuaternionRotation(Q(N_w-1,1:4),[1;0;0]);
D2m0 = QuaternionRotation(Q(N_w-1,1:4),[0;1;0]);
D3m0 = QuaternionRotation(Q(N_w-1,1:4),[0;0;1]);

D1m = QuaternionRotation(Q(N_w,1:4),[1;0;0]);
D2m = QuaternionRotation(Q(N_w,1:4),[0;1;0]);
D3m = QuaternionRotation(Q(N_w,1:4),[0;0;1]);

% Diagonal block

constraint_part = -0.5*DL*rcross(Lambda(:,end))*rcross(D1m);
% constraint_part = 0.5*DL*rcross(Lambda(:,end))*rcross(D1m); %

beta = -strain_twist(3,end) + 0.5*(D2m' * D3m0 - D2m0' * D3m)/DL;

elastic_part = -0.5*(D1m0 + D1m)*(cross(D3m0,D2m) - cross(D2m0,D3m))'/DL;

elastic_part = -0.5*KT*(beta*rcross(D1m) + elastic_part);

elastic_part = elastic_part - KB*(rcross(D1m0)'*rcross(D1m)/DL...
    - 0.5*strain_twist(1,end)*rcross(D2m)...
    - 0.5*strain_twist(2,end)*rcross(D3m));

dTm_dun = Tfac*(elastic_part + constraint_part);

Lie_mat = rcross(U(1:3,N_w))';
Omega_mat = rcross(Omega(:,N_w));

block = -2*dt*(dTm_dun - 0.5*(Omega_mat + Lie_mat*dTm_dun) + ...
    (rcross(cross(U(1:3,N_w),Omega(:,N_w))) + Lie_mat*Omega_mat + Lie_mat^2*dTm_dun)/12)/3;

E(end-2:end,end-2:end) = block;

% Off-diagonal block

elastic_part = 0.5*(D1m0 + D1m)*(cross(D3m0,D2m) - cross(D2m0,D3m))'/DL;

elastic_part = -0.5*KT*(beta*rcross(D1m0) + elastic_part);

elastic_part = elastic_part - KB*(rcross(D1m)*rcross(D1m0)/DL...
    - 0.5*strain_twist(1,end)*rcross(D2m0)...
    - 0.5*strain_twist(2,end)*rcross(D3m0));

dTm_dun = Tfac*elastic_part;

block = -2*dt*(dTm_dun - 0.5*Lie_mat*dTm_dun + Lie_mat^2*dTm_dun/12)/3;

E(end-2:end,end-5:end-3) = block; % Placed at n=m-1

% Now for the interior points, which experience an elastic moment and a
% constraint torque on each side, giving their entries the most general
% form.

for m=2:N_w-1
    
    Tfac = 1/(8*pi*mu*R(m)^3);
    
    D1m0 = QuaternionRotation(Q(m-1,1:4),[1;0;0]);
    D2m0 = QuaternionRotation(Q(m-1,1:4),[0;1;0]);
    D3m0 = QuaternionRotation(Q(m-1,1:4),[0;0;1]);
    
    D1m = QuaternionRotation(Q(m,1:4),[1;0;0]);
    D2m = QuaternionRotation(Q(m,1:4),[0;1;0]);
    D3m = QuaternionRotation(Q(m,1:4),[0;0;1]);
    
    D1m1 = QuaternionRotation(Q(m+1,1:4),[1;0;0]);
    D2m1 = QuaternionRotation(Q(m+1,1:4),[0;1;0]);
    D3m1 = QuaternionRotation(Q(m+1,1:4),[0;0;1]);
    
    Lie_mat = rcross(U(1:3,m))';
    Omega_mat = rcross(Omega(:,m));
    
    % Construct the diagonal block first.
    
    constraint_part = -0.5*DL*rcross(Lambda(:,m-1) + Lambda(:,m))*rcross(D1m);
    % constraint_part = 0.5*DL*rcross(Lambda(:,m-1) +
    % Lambda(:,m))*rcross(D1m); %
    
    beta = -strain_twist(3,m) + 0.5*(D2m1' * D3m...
        - D2m' * D3m1)/DL;
    
    elastic_part_right = 0.5*(D1m + D1m1)...
        *(cross(D3m,D2m1)...
        - cross(D2m,D3m1))'/DL;
 
    elastic_part_right = 0.5*KT*(beta*rcross(D1m) + elastic_part_right);
    
    elastic_part_right = elastic_part_right + KB*(rcross(D1m1)*rcross(D1m)/DL...
        - 0.5*strain_twist(1,m)*rcross(D2m)...
        - 0.5*strain_twist(2,m)*rcross(D3m));

    beta = -strain_twist(3,m) + 0.5*(D2m' * D3m0...
        - D2m0' * D3m)/DL;
    
    elastic_part_left = -0.5*(D1m0 + D1m)...
        *(cross(D3m0,D2m)...
        - cross(D2m0,D3m))'/DL;

    elastic_part_left = -0.5*KT*(beta*rcross(D1m) + elastic_part_left);
    
    elastic_part_left = elastic_part_left - KB*(rcross(D1m0)'*rcross(D1m)/DL...
        - 0.5*strain_twist(1,m)*rcross(D2m)...
        - 0.5*strain_twist(2,m)*rcross(D3m));
    
    dTm_dun = Tfac*(constraint_part + elastic_part_right + elastic_part_left);
    
    block = -2*dt*(dTm_dun - 0.5*(Omega_mat + Lie_mat*dTm_dun) + ...
        (rcross(cross(U(1:3,m),Omega(:,m))) + Lie_mat*Omega_mat + Lie_mat^2*dTm_dun)/12)/3;
    
    E(3*(m-1)+1:3*(m-1)+3,3*(m-1)+1:3*(m-1)+3) = block;
    
    % Now for the off-diagonal blocks, starting with the left one.
    
    elastic_part = 0.5*(D1m0 + D1m)...
        *(cross(D3m0,D2m)...
        - cross(D2m0,D3m))'/DL;
   
    elastic_part = -0.5*KT*(beta*rcross(D1m0) + elastic_part);
    
    elastic_part = elastic_part - KB*(rcross(D1m)*rcross(D1m0)/DL...
        - 0.5*strain_twist(1,m)*rcross(D2m0)...
        - 0.5*strain_twist(2,m)*rcross(D3m0));

    dTm_dun = Tfac*elastic_part;
    
    block = -2*dt*(dTm_dun - 0.5*Lie_mat*dTm_dun + Lie_mat^2*dTm_dun/12)/3;
    
    E(3*(m-1)+1:3*(m-1)+3,3*(m-2)+1:3*(m-2)+3) = block; % Placed at n=m-1
    
    % Finally, the right off-diagonal block.
    
    elastic_part = -0.5*(D1m + D1m1)...
        *(cross(D3m,D2m1)...
        - cross(D2m,D3m1))'/DL;

    elastic_part = 0.5*KT*(beta*rcross(D1m1) + elastic_part);
    
    elastic_part = elastic_part + KB*(rcross(D1m)'*rcross(D1m1)/DL...
        - 0.5*strain_twist(1,m)*rcross(D2m1)...
        - 0.5*strain_twist(2,m)*rcross(D3m1));
    
    dTm_dun = Tfac*elastic_part;
    
    block = -2*dt*(dTm_dun - 0.5*Lie_mat*dTm_dun + Lie_mat^2*dTm_dun/12)/3;
    
    E(3*(m-1)+1:3*(m-1)+3,3*(m)+1:3*(m)+3) = block; % Placed at n=m+1
    
end

% All that remains is to add on the identity part.

for i=1:3*N_w
    
    E(i,i) = E(i,i) + 1;
    
end

end
