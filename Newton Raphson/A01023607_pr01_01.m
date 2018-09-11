% Alberto Pascal A01023607
% Saúl E. Labra A01020725
% Rodrigo García A01024595
% Manuel Guadarrama A01020829

clc;
clear all;
close all;

%Storing the constants of the problem
R=0.08205;
T=273.15;
pressures=[1,2,5,10,40,60,80,100];
V=(R*T)./pressures;
P=0;
Bo=0.05587;
A0=2.2789;
c=128000;
a=0.01855;
b=-0.01587;
% Other variables;
B= (R.*T.*Bo)-A0-((R.*c)./(T.^2));
G= ((-1.*(R.*T.*Bo.*b)) + (A0.*a) -((R.*c.*Bo)./(T.^2)));
D= (R.*Bo.*b.*c)./(T.^2);

num_iter=1;
Initial_Vol=zeros(1,2*length(pressures));
%Arrays to store results
Vol_end=zeros(1,2*length(pressures));
Err_End=zeros(1,2*length(pressures));

%Defining the function and its derivative
 f= @(V, x)((R*T)./V) + (B./(V.^2)) + (G./(V.^3))+(D./(V.^4)) -pressures(x);
 %Array to store results of the function
 Arr_fun=zeros(16,20);
 df=@(V) -1.*((2.*B.*(V.^2) + (4.*D) + (3.*G.*V) + (R.*T.*(V.^3)))./(V.^5));

%Defining tolerance and intial error
Err_tol=.000001;
error = 1000000;
%Array to store Z (compresibility factor)
z_zero=zeros(1,length(pressures));
z_two=zeros(1,length(pressures));

%Iterates over each scenario (pressure)
for x=1:1:length(pressures)
   P=pressures(x); 
   Volum= R*T/P;
   Initial_Vol(1,x)=Volum;
   %Newton Raphson is executed while the iterations doesn't surpass 20 and
   %the error tolerance is fulfilled
   while (error>Err_tol && num_iter<20)
      
       %The function is evaluated and saved with the current pressure iteration
       num1=f(Volum, x);
       Arr_fun(x, num_iter)=num1;
       %The derivative of the function is evaluated with the current pressure iteration
       num2=df(Volum);
       %Newton Raphson formula is applied
       Vnew = Volum - (num1/num2);
       %After the first iteration error can be calculated
       if num_iter>1
           error =(100*(Vnew-Volum)/Vnew);
       end
       Volum=Vnew;
       num_iter=num_iter+1;
   end
   %Volume and error after NR is saved
   Vol_end(x)=Volum;
   Err_End(x)=error;
   %Compresibility factor is calculated
   z_zero(x) = P*Volum/(R*T);
   disp(['The volume ' num2str(x) ' is ' num2str(Vol_end(x)) ' with an error of ' num2str(Err_End(x)) ' after iteration ' num2str(num_iter -1)]);

    %Resetting initial error
   error=1000000;
   num_iter=1;
end

%The same code is executed for the second temperature
T2=473.15;
B= (R.*T2.*Bo)-A0-((R.*c)./(T2.^2));
G= ((-1.*(R.*T2.*Bo.*b)) + (A0.*a) -((R.*c.*Bo)./(T2.^2)));
D= (R.*Bo.*b.*c)./(T2.^2);
f2= @(V, x)((R*T2)./V) + (B./(V.^2)) + (G./(V.^3))+(D./(V.^4)) -pressures(x);
df2=@(V) -1.*((2.*B.*(V.^2) + (4.*D) + (3.*G.*V) + (R.*T2.*(V.^3)))./(V.^5));

for x=1:1:length(pressures)
   P=pressures(x); 
   Volum= R*T2/P;
   Initial_Vol(1,x+8)=Volum;
   %Newton Raphson
   while (error>Err_tol && num_iter<20)
   
      
       num1=f2(Volum, x);
       Arr_fun(x, num_iter)=num1;
       num2=df2(Volum);
       Vnew = Volum - (num1/num2);
       if num_iter>1
           error =(100*(Vnew-Volum)/Vnew);
       end
       Volum=Vnew;
       num_iter=num_iter+1;
   end
   Vol_end(x+8)=Volum;
   z_two(x) = P*Volum/(R*T2);
   Err_End(x+8)=error;
   error=1000000;
   disp(['The volume ' num2str(x+8) ' is ' num2str(Vol_end(x+8)) ' with an error of ' num2str(Err_End(x+8)) ' after iteration ' num2str(num_iter -1)]);
   num_iter=1;
end

%A plot of the compresibility factors vs the pressure is shown
plot(pressures,z_zero,'m^', pressures,z_two, 'r*');
title('Pressure vs Compression values');
ylabel('Compression values [NU]') % x-axis label
xlabel('Pressure values [atm]')


