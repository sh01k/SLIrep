%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama 2012.06.29
% Differential of spherical Bessel function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = sph_besselj_diff(N,X)

Z = (N.*sph_besselj(N-1,X)-(N+1).*sph_besselj(N+1,X))./(2.*N+1);

end