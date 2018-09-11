%Saúl Enrique Labra Cruz A01020725
%Rodrigo García Mayo A01024595
%Alberto Pascal Garza A01023607
%Manuel Guadarrama Villarroel A01020829

clc;
clear all;
close all;

k = 1500;
b = 30;
m = 10;
F_x = -100;
g = -9.81;
F_0 = m * g;

x0 = 0;
y0 = 0;
z0 = 0;
Fz = @(x, y, z) (F_0 - k*y - b*z)/m;

Fz2=@(x,y,z) (F_0 + F_x - k*y - b*z)/m;

Estado_0 = RK4(x0, y0, z0, Fz);
Estado_1= RK4(Estado_0(end,1),Estado_0(end,2), (Estado_0(end,3)), Fz2);
Estado_2= RK4(Estado_1(end,1),Estado_1(end,2), (Estado_1(end,3)), Fz);
Estado_2(end,:)=[];
%Estado_F = RK4(arr(size(arr,1),1),arr(size(arr,1),2)
disp(['Para estabilizarse toma:']);
disp(['Al poner la masa: ' num2str(Estado_0(end,1)) ' segundos']);
disp(['De el punto anterior al aplicar la fuerza: ' num2str(Estado_1(end,1)-Estado_0(end,1)) ' segundos']);
disp(['De el punto anterior a soltar la fuerza: ' num2str(Estado_2(end,1)-Estado_1(end,1)) ' segundos']);

disp(['Total: ' num2str(Estado_2(end,1)) 'segundos']);


plot(Estado_0(:,1), Estado_0(:,2)), hold on
plot(Estado_1(:,1), Estado_1(:,2)), hold on
plot(Estado_2(:,1), Estado_2(:,2));
title("Time vs Position");
xlabel("Time [s]");
ylabel("Position [m]");
legend("With mass only", "Added Force", "Released state");
figure;
plot(Estado_0(:,1), Estado_0(:,3)), hold on
plot(Estado_1(:,1), Estado_1(:,3)), hold on
plot(Estado_2(:,1), Estado_2(:,3));
title("Time vs Velocity");
xlabel("Time [s]");
ylabel("Velocity [m/s]");
legend("With mass only", "Added Force", "Released state");


figure;

plot(Estado_0(:,1), Estado_0(:,4)), hold on
plot(Estado_1(:,1), Estado_1(:,4)), hold on
plot(Estado_2(:,1), Estado_2(:,4));
title("Time vs Acceleration");
xlabel("Time [s]");
ylabel("Acceleration [m/s^2]");
legend("With mass only", "Added Force", "Released state");



%1=TIEMPO, 2=POSICIÓN 3=VELOCIDAD 4=ACELERACIÓN 5=ENERGÍA POTENCIAL 6=
%ENERGÍA  7= total energy;
maxVE0 = max(Estado_0(:,3)); % - min(Estado_0(:,3));
minVE0 = min(Estado_0(:,3)); % - min(Estado_2(:,3));
maxVE1 = max(Estado_1(:,3)); % - min(Estado_0(:,3));
minVE1 = min(Estado_1(:,3)); % - min(Estado_2(:,3));
maxVE2 = max(Estado_2(:,3)); % - min(Estado_0(:,3));
minVE2 = min(Estado_2(:,3)); % - min(Estado_2(:,3));

Estado_2=Energy_Calculations(Estado_2, m, g,k, Estado_1(end,2));
Estado_1 = Energy_Calculations(Estado_1, m, g, k, Estado_0 (end,2));

Lost_between_1_2 = Lost_Energy(Estado_1(end,7), Estado_2);


figure;
plot(Estado_2(:,1),Estado_2(:,5)), hold on;
plot(Estado_2(:,1), Estado_2(:,9)), hold on;
plot(Estado_2(:,1), Estado_2(:,6));
title("Energies during the release state");
xlabel("Time [s]");
ylabel("Energy [J]");
legend("Potential-Gravitational", "Potential-Elastic", "Kinetic");


%figure;
%plot(Estado_2(:,1),Estado_2(:,6));
%title("Kinetic Energy during the release state vs Time");
figure;
plot(Estado_2(:,1), Estado_2(:,7));
title ("Total Energy during the release state (Sum of the energies)");
xlabel("Time [s]");
ylabel("Energy [J]");
figure;
plot(Estado_2(:,1), Lost_between_1_2(:,1));
title("Lost Energy");
xlabel("Time [s]");
ylabel("Energy [J]");
%figure;
%plot(Estado_1(:,1), Estado_1(:,7));
%title("Total energy from state 1");

[a,b] = Regression(Estado_2(:,1), Estado_2(:,7), 33600);

func = @(x) a*exp(b*x);
figure; 
plot(func(0:20)), hold on; 
title("Regression of the damping");
xlabel("Time [s]");
ylabel("Energy [J]");
func2 = @(x) a*exp(b*x)+Estado_2(end,7);
plot(func2(0:20)); 
legend("Original Regression (23.6019 e^-^0^.^3^9^9^4^x)", "Added last energy value from previous state(23.6019 e^-^0^.^3^9^9^4^x) + [9.7933]");
%We are not completely sure about these damping units above

[Estado_0, Estado_1, Estado_2] = Displacement (Estado_0, Estado_1, Estado_2);
figure;
plot(Estado_0(:,1), Estado_0(:,2)), hold on
plot(Estado_1(:,1), Estado_1(:,2)), hold on
plot(Estado_2(:,1), Estado_2(:,2)), hold on
title("Time vs Position (Displaced) With Roots as O");
xlabel("Time [s]");
ylabel("Position (Displaced) [m]");
%legend("Mass Only", "Added Force", "Released State", "Roots");





%Calculating the crossing }

Roots = Find_Roots(Estado_2);

function [root] = Find_Roots(Estado_2)
    root=[];
    r=1;
    index = 0;
    p = 0;
    t = 0;
    while(index+1 >0 && index+1< size(Estado_2,1))
        num = Newton_Raphson (t, p, Estado_2(index+1,3));
        root = [root;num];
        plot(num, 0, "o"), hold on;
        [t, p, index] = Sign (Estado_2, index+1);
        r=r+1;
    end
    root(1)=[];
    
end
function [xi] = Newton_Raphson (xi, fxi, dfxi)
error = 100;
n = 0;
while (error < 0.00000000001)
    xnew  = xi - (fxi/dfxi);
    n = n+1;
    if(n>1)
        error= abs(((xnew-xi)/xnew));
    end
    xi = xnew;
end
end

function [t, p, index] = Sign (Estado_2, posi)
xnext = 0;
p = 0;
t=0;
index = posi;
for i=posi:size(Estado_2,1)-1
        xnext = Estado_2(i+1,2);
        if (Estado_2(i,2) < 0 && xnext >0)
            p = xnext;
            t = Estado_2(i,1);
            index = i;
            break;
        elseif (Estado_2(i,2) > 0 && xnext <0)
            p = xnext;
            t = Estado_2(i,1);
            index = i;
            break;
        else
            index=-100;
        end     
end
end
function [Estado_0, Estado_1, Estado_2] = Displacement (Estado_0, Estado_1, Estado_2)
    move = Estado_0(size(Estado_0, 1),2);
    for i=1:size(Estado_0,1)
        Estado_0(i,2) = Estado_0(i,2) - move;
        Estado_1(i,2) = Estado_1(i,2) - move;
        Estado_2(i,2) = Estado_2(i,2) - move;
    end
end
function [Lost] = Lost_Energy(Energy, Estado)
for n=1: size(Estado(:,7))
    Lost(n,1) = Energy - Estado(n,7);
end
end
function [Estado] = Energy_Calculations(Estado, mass, gravity, k, previous_position)
for n=1:size(Estado(:,2))
   
%Potential Energy 
%(Gravitational) + (Elastic)
Estado(n,5) = mass*gravity*Estado(n,2);
%Elastic Energy
Estado(n,9)=((0.5)*k*((Estado(n,2)-previous_position)^2));
%Kinetic Energy
Estado(n,6) = 0.5*mass*(Estado(n,3)^2);

Estado(n,7) = Estado(n,6) + Estado (n,5) + Estado(n,9);

end


end
function [arr] = RK4(x0, y0, z0, Fz) 
    xactual = x0;
    yactual = y0;
    zactual = z0;
    avg_y = 1000000;
    arr = [];
    cont = 1;
    dx = 0.0001;
    fx1 = @(y,k1,k2,k3,k4) (y + (1/6)*(k1+2*k2+2*k3+k4)*dx);
    while (avg_y > 1e-5)
        k1_y = zactual;
        k1_z = Fz (xactual, yactual, zactual);
        k2_y = zactual + 0.5 * k1_z * dx;
        k2_z = Fz (xactual + 0.5*dx, yactual+ 0.5 * k1_y * dx, zactual + 0.5 * k1_z * dx);
        k3_y  = zactual + 0.5 * k2_z * dx;
        k3_z = Fz (xactual + 0.5*dx, yactual+ 0.5 * k2_y * dx, zactual + 0.5 * k2_z * dx);
        k4_y = zactual + k3_z * dx;
        k4_z = Fz (xactual + dx, yactual+ k3_y * dx, zactual + k3_z * dx);
        yactual = fx1 (yactual, k1_y, k2_y, k3_y, k4_y);
        xactual = xactual + dx;
        zactual = fx1 (zactual, k1_z, k2_z, k3_z, k4_z);
        arr(cont,1) = xactual;
        arr(cont, 2) = yactual;
        arr(cont, 3) = zactual;
        arr(cont, 4) = Fz(xactual, yactual, zactual);
        if (cont > 5)
            avg_y = mean(abs(arr(cont-4:cont,3)));           
        end
        cont = cont + 1;
    end
    %disp(cont);
end
function [A0, A1] = Regression (Times, Energies, n)
A0=1;
A1=1;
%Calculating the necesary changes:
Ln_Times=[];
Ln_Energies=[];

for x=1:n
    Ln_Times(x,1) = log(Times(x,1));
    Ln_Energies(x,1) = log(Energies(x,1));
end

%Calculating the sumatories
Sum_LnTimes=0;
Sum_LnEnergies=0;
Sum_LnTimesxEnergies=0;
Sum_LnSquareTimes=0;
for t=1:n
   Mult(t,1)=Ln_Times(t,1) * Ln_Energies(t,1); 
end
for i=1:n
   Sum_LnTimes= Sum_LnTimes + Ln_Times(i,1);
   Sum_LnEnergies= Sum_LnEnergies + Ln_Energies(i,1);
   Sum_LnTimesxEnergies = Sum_LnTimesxEnergies + (Ln_Times(i,1) * Ln_Energies(i,1));
   Sum_LnSquareTimes = Sum_LnSquareTimes + (Ln_Times(i,1)^2);
end

Q1= sum(Ln_Times,1);
Q2=sum(Ln_Energies,1);
Q3=sum(Mult,1);
%we prepare the matrix:
MatA = [n, Sum_LnTimes;Sum_LnTimes, Sum_LnSquareTimes];
MatB = [Sum_LnEnergies; Sum_LnTimesxEnergies];

Sols=MatA\MatB;

%Now we convert back
A0 = exp(Sols(1,1));
A1=Sols(2,1);
end