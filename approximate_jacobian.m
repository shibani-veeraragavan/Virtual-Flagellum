function [J0] = approximate_jacobian(Filament,dt,mu)
% APPROXIMATE_JACOBIAN  Creates the Jacobian for a given filament based on
%                       analytic expressions given in Appendix B of the
%                       paper.

N_w = Filament.N_w; % The number of particles in the filament.

J0 = zeros(6*N_w);

% We start with the entries of the Jacobian corresponding to the update
% equation for the position of the first particle.
j = 3*N_w + 3;
fac = dt/(9*pi*mu*Filament.R(1));

J0(1,1) = 1; J0(2,2) = 1; J0(3,3) = 1;
J0(1,j+1) = fac; J0(2,j+2) = fac; J0(3,j+3) = fac;

% Add [I] row beside [D] %
% for I=2:N_w
%     m=(I-1)*3;
%     J0(m+1:m+3,1:3) = eye(3);
% end
% Next, we provide the terms related to the derivative of the velocity
% constraints with respect to the Lie algebra elements.
% Why divide by -2dt/3 ? Meant to be divided to mu at J_block_D.
J0(4:3*N_w,4:3*N_w+3) = J_block_C(Filament,-0.75*Filament.DL/dt);
% J0(4:3*N_w,4:3*N_w+3) = J_block_C(Filament,0.5*Filament.DL); %

% Now, we produce the block relating the velocity constraints to the
% Lagrange multipliers.

J0(4:3*N_w,3*N_w+4:end) = J_block_D(Filament,dt,mu);

% Next, the block encoding the dependence of the Lie algebra update
% equations on the Lie algebra elements.

J0(3*N_w+1:end,4:3*N_w+3) = J_block_E(Filament,dt,mu);

% Finally, we produce the derivatives of the Lie algebra update equations
% with respect to the Lagrange multipliers.

J0(3*N_w+1:end,3*N_w+4:end) = J_block_F(Filament,dt,mu);

end