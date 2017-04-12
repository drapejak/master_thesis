function G=gaussMean(x,sigma,mean)
G = exp(-(x - mean).^2/(2*sigma^2))/(sqrt(2*pi)*sigma);