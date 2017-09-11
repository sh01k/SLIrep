%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama 2012.06.29
% Differential of spherical Hankel function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = sph_besselh_diff(N,K,X)

Z = (N.*sph_besselh(N-1,K,X)-(N+1).*sph_besselh(N+1,K,X))./(2.*N+1);

end