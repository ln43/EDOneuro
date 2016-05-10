function y=morrislecar(t,y,par)
Minf=1/2*(1+tanh((y(1)-par(9))/par(10)));
Ninf=1/2*(1+tanh((y(1)-par(11))/par(12)));
lambN=1/par(14)*cosh((y(1)-par(11))/(2*par(12)));

y=[1/par(13)*(par(1)-par(3)*(y(1)-par(6))-par(4)*Minf*(y(1)-par(7))...
    -par(5)*y(2)*(y(1)-par(8)));lambN*(Ninf-y(2))];