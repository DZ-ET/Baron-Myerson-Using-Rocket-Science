function x = GG(y,par)
    x = interp1(par.v, par.G, y, 'spline', 'extrap');
    x(x>1) = 1;
    x(x<0) = 0;
end
