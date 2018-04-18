function [fx] = pdfAtX(pdf, mu, sig, x)
% Returns the probability density "fx" of the values in "x" (scalar, vector
% or matrix), according to the probability density "pdf", specified by mean
% "mu" and standard deviation "sig".
%
% Syntax:
% [fx] = pdfAtX(pdf, mu, sig, x)
%
% pdf: string for probability density function among... 
% %......n, Gaussian
% %......w, Inverse Gaussian (Wald)
% %......l, Lognormal
% %......g, Gamma
% %......i, Inverse-Gamma/Scaled-inverse-chi-square
% %......e, Exponential... only valid when the mean is equal to the
%           standard deviation!

switch pdf
    case 'n'% Normal
        fx = (1/sqrt(2*pi*sig^2))*exp(-((x - mu).*(x - mu))/(2*sig^2));
        
    case 'w'% Inverse Gaussian
        lambda=mu^3/sig^2;
        fx=sqrt(lambda./(2*pi*x.^3)).*...
            exp((-lambda*(x-mu).^2)./(2*(mu^2)*x));
        
    case 'l'% Lognormal
        lognMu = log(mu^2/sqrt(sig^2 + mu^2));
        lognSig =  sqrt(log(((sig/mu)^2) + 1));
        fx = (1./(x*lognSig*sqrt(2*pi))).*...
            exp(-((log(x) - lognMu).*(log(x) - lognMu))/(2*lognSig^2));
        
    case 'g'% Gamma
        alpha = (mu/sig)^2;
        beta = mu/(sig^2);
        fx = ((beta^alpha*x.^(alpha-1))/gamma(alpha)).*exp(-beta*x);
        
    case 'i'% Inverse-Gamma/Scaled-Inverse-Chi-Square
        alpha = ((mu/sig)^2) + 2;
        beta = (mu^3/sig^2) + mu;
        fx=((beta^alpha*x.^(-alpha-1))/gamma(alpha)).*exp(-beta./x);
        
    case 'e'% Exponential, using only the mean!
        lambdae = 1/mu;
        fx = lambdae*exp(-lambdae.*x);
        
    otherwise % Lognormal densities by default...
        lognMu = log(mu^2/sqrt(sig^2 + mu^2));
        lognSig =  sqrt(log(((sig/mu)^2) + 1));
        fx = (1./(x*lognSig*sqrt(2*pi))).*...
            exp(-((log(x) - lognMu).*(log(x) - lognMu))/(2*lognSig^2));
        warning(['Density "' pdf '" not recognized,'...
            ' lognormal densities returned by default']);
end