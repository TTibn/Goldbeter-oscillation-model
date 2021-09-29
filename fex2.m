## Copyright (C) 2020 Theocharis
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} f_test (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Theocharis <Theocharis@THEOCHARIS-HP>
## Created: 2020-11-15

function pdot = fex2(p,t)
    global vs;
    global vm;
    global vd;
    global ks;
    global k1;
    global k2;
    global V1;
    global V2;
    global V3;
    global V4;
    global K1;
    global K2;
    global K3;
    global K4;
    global KI;
    global Km1;
    global Kd;
    global n;
    pdot=zeros(5,1);
    pdot(1)=vs*((KI^n)./((KI^n)+((p(5)^n))))- vm*((p(1))./(Km1+p(1)));                     
    pdot(2)=(ks*p(1))-V1*(p(2)./(K1+p(2)))+(V2*(p(3)./(K2+p(3))));
    pdot(3)=V1*(p(2)./(K1+p(2)))-(V2*(p(3)./(K2+p(3))))-(V3*(p(3)./(K3+p(3))))+(V4*(p(4)./(K4+p(4))));
    pdot(4)=V3*(p(3)./(K3+p(3)))-(V4*(p(4)./(K4+p(4))))-(k1*p(4))+(k2*p(5))-(vd*(p(4)./(Kd+p(4))));
    pdot(5)=k1*p(4)-(k2*p(5));
endfunction
