% @fileName: randInvGauss.m
% @purpose: generate random numbers distributed as an inverse-Gaussian
%                specified by the paramters in the input.
%
% @syntax:
%           [x] = randInvGauss(mu, lambda, noObs, noCh)
% 
% where: 
% mu: mean
% lambda: natural shape parameter of the inverse Gaussian (mean^3/std^2)
% noObs: number of observations (rows for the matrix returned)
% noCh: number of parallel series (columns)
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

function [x] = randInvGauss(mu, lambda, noObs, noCh)
% random Gaussian data 
normRandData = random('norm', 0, 1, noObs, noCh);
y = normRandData.*normRandData;
x = mu + ((mu^2 * y)/(2*lambda)) -...
    ((mu/(2*lambda))*(sqrt((4*mu*lambda*y)+(mu^2*y.*y))));
z = random('unif', 0, 1, noObs, noCh);
for row = 1:noObs
    for column = 1:noCh
        if z(row, column) > mu./(mu+x(row, column))
            x(row,column) = (mu^2)/x(row, column);
        end
    end
end
