function [ g, err, b ] = ShrinkwrapFunction( diff_pattern, n )
% SHRINKWRAPFUNCTION Shrinkwrap algorithm
% Susanna Hammarberg 20170129
% Following Reference [1] Marchesini et al.,
% "X-ray image reconstruction from a diffraction pattern alone"
% Phys. Rev. B 68, 140101 (2003)
% Function needed to runt this code: create2Dgaussian.m
% Input:
% diff_patter - diffraction pattern as 2D matrix
% n - nbr of iterations
% Output:
% g - reconstructed image
% err - error evaluation at every iteration
% b - the best iteration run.

% the intensity of the diffraction pattern is F^2
absF = diff_pattern.^0.5;

% parameters
beta = 0.9;              % HIO feedback parameter
thresholdMask1 = 0.04;   % Initial threshold for mask
thresholdMask2 = 0.2;    % threshold for mask 
sigma = 3;               % initial value of gaussian width
c = size(absF,1)/2;      % center of Gaussian

% initialize vector for error calculation
err2 = zeros(1,n);

% g represents the unknown object. Create initial guess for g using 
% the autocorrelation function of the object. Start by Fourier 
% transform the diffraction pattern:
g = fftshift(ifft2(absF,'symmetric'));
% create autocorrelation function of g
g = abs(fftshift(ifft2(fft2(g).*conj(fft2(g)))))./(size(g,1)*size(g,2));
% create a logical matrix that is the Shrinkwrap support
support = g > thresholdMask1*max(g(:));
% define iteration counter
k = 0;  

% start iterations
while k < n
    % every 20th iteration, update the support:
    if mod(k,20)==0 && k>0
        % call function create2Dgaussian
        gauss2D = create2Dgaussian(sigma,c,size(absF,1),size(absF,2));
        % calculate the convolution the absolute value of the object wave 
        % with a gaussian
        support = conv2(gauss2D,abs(g),'same');
        % Create a logical matrix that is the Shrinkwrap support
        support = support > thresholdMask2.*max(support(:));        
        % reduce sigma with 1% every 20th iteration until a limit of 1.5
        if sigma >= 1.5
            sigma = sigma*0.99;
        end
%         uncomment to display support at every support update
%         figure
%         imagesc(support)
    end
        
    % STEP 1: Fourier transform of g(x),  G(u) :
    G = fftshift(fft2(g));  
    
    % STEP 2: make |G| confirm with |F|
    Gprime = absF.*exp(1i.*angle(G)); 
    
    % STEP 3: inverse Fourier transform      
    gprime = (ifft2(ifftshift(Gprime)));
   
    % STEP 4: See to that gprime satisfies its constraints: 
     for l= 1:size(absF,1)
        for j=1:size(absF,2) 
            if support(l,j) == 0   
                % update g(l,j) outside the support
                gprime(l,j) = g(l,j) - beta .* gprime(l,j);
                
            else
                % update g'(l,j) inside the support in one of 
                % the two following ways:
                %gprime(l,j) = g(l,j);
                gprime(l,j) = gprime(l,j);
            end
        end
     end

    % overwrite g with result from iteration
    g = gprime;
    % set all negative values of g to zero
    % NOT following [1]
    g(g<0) = 0;  
    
    k = k+1;
   
%     %uncomment to display image at every iteration:   
%     figure(k)    
%     imshow((g))
%     hold on
%     title(k)
    
    % calculate error of iteration
    err2(k) = sum(sum( (abs(abs(fftshift(fft2(((g.*conj(g))))))) - absF).^2))...
        ./ sum(sum(absF.^2));
    
end %End of iterations

% Figure plotting
figure(9987)
err = err2.^0.5;
plot(2:k,err(2:end))
title('Error per iteration')
xlabel('k')
ylabel('Error')
% best result at itertion 'b'
[a b]=min(err);
disp('The best result is at iteration:')
disp(b)

end

