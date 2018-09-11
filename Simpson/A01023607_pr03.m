%Saúl Enrique Labra Cruz A01020725
%Rodrigo García Mayo A01024595
%Alberto Pascal Garza A01023607
%Manuel Guadarrama Villarroel A01020829

clear *;
close all;
clc;

%Load data from external files
load("batrang.mat");
x = xx;
y = yy;
load("batrangBroken.mat");
x_b = xx;
y_b = yy;

dx = 0.01;

%Displace the graph upwards, in order to have all points above x-axis
ynewup = Desplazamiento(y);
y_bnewup = Desplazamiento(y_b);
xnew = Desplazamiento(x);

plot (x_b,y_bnewup,'o'), hold on
xlabel('Units: centimeters');
ylabel('Units: centimeters');
plot (x,ynewup,'-'), hold on

%Saves the values into two separate arrays representing the two curves that
%conform the shape
[ynewdown, ynewup] = F_Div(x,ynewup);
[y_bnewdown, y_bnewup] = F_Div(x_b,y_bnewup);

% plot(x,y,'*-');
% plot(x(1:1395),ynewup,'*-');
% plot(x_b,y_b,'*-');
% plot(x_b,y_bnewup,'*-');

%Gets the integral of the 4 curves
Integral_up = Simpson (ynewup,dx);
Integral_b_up = Simpson (y_bnewup,dx);
Integral_down = Simpson (ynewdown,dx);
Integral_b_down = Simpson (y_bnewdown,dx);

%Gets the area of the two shapes (batrang original, batrang broken)
Integral = double(Integral_up - Integral_down);
Integral_b = double(Integral_b_up - Integral_b_down);

%Gets the center of mass of the broken batrang

%Center in X: Applies the formula multiplication to get the center of mass
%before integrating
fx_tosend= x_b(1:length(y_bnewup)).*y_bnewup;
gx_tosend = x_b(length(y_bnewup)+1:end).*y_bnewdown;
CentroX_b = Simpson_centrosX (fx_tosend, gx_tosend, dx);
CentroX_b = double(CentroX_b / Integral_b);

%Center in Y: Applies the formula "2" exponent to get the center of mass
%before integrating
fx_tosend= y_bnewup.^2;
gx_tosend= y_bnewdown.^2;
CentroY_b = Simpson_centrosY (fx_tosend, gx_tosend, dx);
CentroY_b = double(CentroY_b / Integral_b);

%Gets the center of mass of the original batrang

%Center in X: Applies the formula multiplication to get the center of mass
%before integrating
fx_tosend= x(1:length(ynewup)).*ynewup;
gx_tosend = x(length(ynewup)+1:end).*ynewdown;
CentroX = Simpson_centrosX (fx_tosend, gx_tosend, dx);
CentroX = double(CentroX / Integral);

%Center in Y: Applies the formula "2" exponent to get the center of mass
%before integrating
fx_tosend= ynewup.^2;
gx_tosend= ynewdown.^2;
CentroY = Simpson_centrosY (fx_tosend, gx_tosend, dx);
CentroY = double(CentroY / Integral);

%Plotting the graph and the values on console

plot(CentroX,CentroY,'*'), hold on
plot(CentroX_b,CentroY_b,'*'), hold on
legend('broken batrang', 'complete batrang', 'complete batrang mass center', 'broken batrang mass center');
title('Center of mass comparison: Batrang');

disp(['Original batrang integral: ' num2str(Integral) ' cm^2 ']);
disp(['Broken batrang integral: ' num2str(Integral_b) ' cm^2']);
disp(['Original batrang mass center: (' num2str( CentroX ) ', ' num2str(CentroY) '). These coordinates are referencing centimeters from the origin']);
disp(['Broken batrang mass center: (' num2str(CentroX_b) ', ' num2str(CentroY_b) '). These coordinates are referencing centimeters from the origin']);

%checking the difference in distance

Check_Diff(CentroX,CentroY, CentroX_b, CentroY_b);

%Function to calculate the integral in Y according to the center of mass formula
function [I] = Simpson_centrosY (fx,gx,dx)
    sumOdd= 0;
    sumEven= 0;

    for i=2:2:length(fx)-1
            sumOdd = sumOdd + 0.5*(fx(1,i) - gx(1,i));
    end
    
    for i=3:2:length(fx)-2
            sumEven = sumEven + 0.5*(fx(1,i) - gx(1,i));
    end
    
    I= (dx/3) * (fx(1,1) + (4*sumOdd) + (2*sumEven) + fx(1,length(fx)));
end

%Function to calculate the integral in X according to the center of mass formula
function [I] = Simpson_centrosX (fx,gx,dx)

    sumOdd= 0;
    sumEven= 0;

    for i=2:2:length(fx)-1
            sumOdd = sumOdd + (fx(1,i) - gx(1,i));
    end
    
    for i=3:2:length(fx)-2
            sumEven = sumEven + (fx(1,i) - gx(1,i));
    end
    
    I= (dx/3) * (fx(1,1) + (4*sumOdd) + (2*sumEven) + fx(1,length(fx)));
end

%Checks if the difference of centers of mass is greater than 1
function Check_Diff(x1,y1,x2,y2)

    difference = (((x2-x1)^2) - ((y2-y1)^2))^(.5);
    
    if(difference >1)
        disp(' ');
        disp(['The difference between the centers of mass is ' num2str(difference)]);
        disp(' ');
        disp('Batman can sue China Corp.');
        disp(' ');
        disp('Batman is not too old to use batrangs');

    else
        disp(' ');
        disp(['The difference between the centers of mass is ' num2str(difference)]);
        disp(' ');
        disp('Batman cannot sue China Corp.');
        disp(' ');
        disp('Batman is too old to use batrangs');

    end
end

%Separates the shape into two functions (up and down)
function [temp1, temp2] = F_Div (x,y)
    temp1 = [];
    temp2 = [];
    cont1 = 1;
    cont2 = 1;
    
    %If the values of the array are returning (change direction in x) the
    %function are saved in another array
    for i=2:length(y)
        if (x(1,i) < x(1,i-1))
            temp1(1,cont1) = y(1,i);
            cont1 = cont1 + 1;
        else
            temp2(1,cont2) = y(1,i);
            cont2 = cont2 + 1;
        end
    end
    temp2 = [3 temp2];
end

%Displaces the y values in the graph upwards so the values become positive
function [y] = Desplazamiento (y)
    minimo = abs(min(y));
    for i=1:length(y)
        y(1,i) = y(1,i) + minimo;
    end
end

%Simpson 1/3 method
function [I] = Simpson (fx,dx)
    sumOdd= 0;
    sumEven= 0;

    for i=2:2:length(fx)-1
            sumOdd = sumOdd + fx(1,i);
    end
    
    for i=3:2:length(fx)-2
            sumEven = sumEven + fx(1,i);
    end
    
    I=(dx/3) * (fx(1,1) + (4*sumOdd) + (2*sumEven) + fx(1,length(fx)));

end