% Creates a 2D gaussian with dimensions [m,n] 
% Input parameters:
% sigma : standars deviation
% c: position of peak center
% m: width of function
% n: length of function
% Susanna Hammarberg 20170129
function [gauss2D] = create2Dgaussian( sigma, c, width, length)
% allocate memory for matrix
gauss2D = zeros(width,length);
% create a 2DGaussian
for m=1:width
    for n=1:length
        gauss2D(m,n) = exp(-(((n-c).^2)+(m-c).^2)./(2.*sigma.^2));
    end
end

end