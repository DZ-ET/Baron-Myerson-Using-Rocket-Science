function x = ff(y,par)
    x = interp1(par.v, par.f, y, 'spline', 'extrap');
    x(y>1) = 0;
    x(y<0) = 0;
end