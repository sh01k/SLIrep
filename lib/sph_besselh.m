%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama 2012.06.29
% Spherical Hankel function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = sph_besselh(N,K,X)

Z = sqrt(pi./(2.*X)).*besselh(N+0.5,K,X);

end