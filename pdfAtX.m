% @fileName: pdfAtX.m
% @purpose: return the probability density "fx" of the values in "x" 
%                (scalar, vector, or matrix), according to the probability 
%                density "pdf", specified by mean "mu" and standard 
%                deviation "sig".
%
% @syntax:
%           [fx] = pdfAtX(pdf, mu, sig, x)
%
% where pdf is a string standing for... 
%............n: Gaussian
%............w: inverse Gaussian (Wald)
%............l:   lognormal
%............g:  gamma
%............i:   inverse-gamma/scaled-inverse-chi-square
%............e:  exponential... only the mean will be used!
%  
% @version: 0.1
% @created:  12-Jan-2016       
% @lastModified: 19-Apr-2018   
% @bugs: no known bugs  
% @author: Javier A. Caballero (jcaballeroo)
% @lastEditor: Javier A. Caballero (jcaballeroo)
% 
% @license: MIT license
%  
%  Copyright (c) 2018, Javier A. Caballero
%  
%  Permission is hereby granted, free of charge, to any person obtaining a copy
%  of this software and associated documentation files (the "Software"), to deal
%  in the Software without restriction, including without limitation the rights
%  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%  copies of the Software, and to permit persons to whom the Software is
%  furnished to do so, subject to the following conditions:
%  
%  - The above copyright notice and this permission notice shall be included in all
%  copies or substantial portions of the Software.
% 
%  - Reference to the original source should be given.
%  
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%  SOFTWARE.


function [fx] = pdfAtX(pdf, mu, sig, x)
switch pdf
    case 'n'% Gaussian
        fx = (1/sqrt(2*pi*sig^2))*exp(-((x - mu).*(x - mu))/(2*sig^2));
        
    case 'w'% inverse Gaussian
        lambda=mu^3/sig^2;
        fx=sqrt(lambda./(2*pi*x.^3)).*...
            exp((-lambda*(x-mu).^2)./(2*(mu^2)*x));
        
    case 'l'% lognormal
        lognMu = log(mu^2/sqrt(sig^2 + mu^2));
        lognSig =  sqrt(log(((sig/mu)^2) + 1));
        fx = (1./(x*lognSig*sqrt(2*pi))).*...
            exp(-((log(x) - lognMu).*(log(x) - lognMu))/(2*lognSig^2));
        
    case 'g'% gamma
        alpha = (mu/sig)^2;
        beta = mu/(sig^2);
        fx = ((beta^alpha*x.^(alpha-1))/gamma(alpha)).*exp(-beta*x);
        
    case 'i'% inverse-Gamma/scaled-inverse-chi-square
        alpha = ((mu/sig)^2) + 2;
        beta = (mu^3/sig^2) + mu;
        fx=((beta^alpha*x.^(-alpha-1))/gamma(alpha)).*exp(-beta./x);
        
    case 'e'% exponential, using only the mean!
        lambdae = 1/mu;
        fx = lambdae*exp(-lambdae.*x);
        
    otherwise % lognormal by default...
        lognMu = log(mu^2/sqrt(sig^2 + mu^2));
        lognSig =  sqrt(log(((sig/mu)^2) + 1));
        fx = (1./(x*lognSig*sqrt(2*pi))).*...
            exp(-((log(x) - lognMu).*(log(x) - lognMu))/(2*lognSig^2));
        warning(['Density "' pdf '" not recognized,'...
            ' lognormal densities returned by default']);
end