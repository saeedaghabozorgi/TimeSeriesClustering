function Fourier_distance()

load leleccum; 
 s = leleccum;
 a=s(1:365);
 b=s(101:465);

%   X= cumsum(randn(100,1));
% y = cumsum(randn(100,1));

figure;
hold off;
a = (a-mean(a))/std(a);         % z-normalization
b = (b-mean(b))/std(b); 
plot(a); 
hold on;
plot(b);
 


fa = fft(a)/sqrt(length(a)); 
fb = fft(b)/sqrt(length(b));

d2=dis_euclidean(a,b);
d3=dis_euclidean(fa,fb);
 
 
%  hold on;
% plot(a_top20, 'r');
%  hold on;
%  plot(b_top20, 'r');
 
fa = fft(a); % Fourier decomposition
fa(10:end-10) = 0; % keep first 5 coefficients (low frequencies)
reconstra = real(ifft(fa)); % reconstruct signal

fb = fft(b); % Fourier decomposition
fb(10:end-10) = 0; % keep first 5 coefficients (low frequencies)
reconstrb = real(ifft(fb)); % reconstruct signal

d4=dis_euclidean(fa/sqrt(length(a)),fb/sqrt(length(a)));

 hold on;
plot(reconstra, 'r');
  hold on;
plot(reconstrb, 'r');
end