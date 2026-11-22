function [K,M,ke] = beamFE(L,E,rho,A,I,ne,mtip,massType)
%BEAMFE  Euler-Bernoulli cantilever finite-element model written from scratch.
%   [K,M,ke] = beamFE(L,E,rho,A,I,ne,mtip,massType)
%   2-node Hermite beam elements, DOFs per node [w theta] (deflection, slope).
%   Node 1 is clamped, so its two DOFs are removed. A point mass mtip is
%   added at the free end. massType = 'consistent' (default) or 'lumped'.
%   Returns reduced global stiffness K and mass M, and the element stiffness ke.
if nargin < 8, massType = 'consistent'; end
le = L/ne;  nd = 2*(ne+1);
ke = E*I/le^3 * [ 12     6*le    -12     6*le;
                  6*le   4*le^2  -6*le   2*le^2;
                 -12    -6*le     12    -6*le;
                  6*le   2*le^2  -6*le   4*le^2];
if strcmpi(massType,'lumped')
    me = rho*A*le/2 * diag([1 0 1 0]);
else
    me = rho*A*le/420 * [ 156     22*le    54     -13*le;
                          22*le   4*le^2   13*le  -3*le^2;
                          54      13*le    156    -22*le;
                         -13*le  -3*le^2  -22*le   4*le^2];
end
K = zeros(nd);  M = zeros(nd);
for e = 1:ne                                  % assembly
    idx = 2*e-1 : 2*e+2;
    K(idx,idx) = K(idx,idx) + ke;
    M(idx,idx) = M(idx,idx) + me;
end
M(end-1,end-1) = M(end-1,end-1) + mtip;       % tip mass on last deflection DOF
K = K(3:end,3:end);  M = M(3:end,3:end);      % clamp node 1
end
