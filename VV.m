function V = VV(z,par)

GrossVal = cumtrapz(par.v(end:-1:1) .* par.g) * par.d; % v represents Q

V = interp1(par.v, GrossVal, z, 'spline', 'extrap');

V(z > 1) = 0;
V(z < 0) = 0;

end