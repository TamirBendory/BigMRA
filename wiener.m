function x = wiener(y,z,sigma)
% we try to estimate x from z = x*y + noise 
% everything is in the Fourier domain

x = (z.*conj(y))./ (abs(y).^2 + sigma^2); 