%Saúl Enrique Labra Cruz A01020725
%Rodrigo García Mayo A01024595
%Alberto Pascal Garza A01023607
%Manuel Guadarrama Villarroel A01020829

clear *;
close all;
clc;
%%Up to here we only recovered the lost data values. 
load("satuSignal.mat");
%plot(tsat,fsat,'*-');

%% Calculating delta x

True_Dx=0.01;
Dx=0;
Dx_num=zeros(1,274);
tsat=round(tsat,4);
tsat_toadd=zeros(1,1);
fsat_toadd=zeros(1,1);
for n=2 : size(tsat,2) %This will calculate all our delta x throughout the function.
   
   Dx(1,n-1) =(tsat(1,n) - tsat(1,n-1)); %our values will be stored in this array. 
   
end
Dx_num = round(Dx,4); %to prevent matlab from doing wierd stuff, we round up to four decimal points. (Values like .0100000 that looked equal were taken as different)
step = mode(Dx_num); %our step (delta x) will be the value that is repeated the most throught the delta xs.


%% Retrieving missing points

%polynomial regression using between1, between2 and the previous
%calculated value.

for n=1: size(Dx_num,2)

    %we will check where our step is bigger than it should.
    if(Dx_num (1, n) > step)
       %if it is bigger, it means we are missing data here. We need to know how many
       missing_points = (Dx_num(1,n)/step ) -1; %%this will tell us how many points we're missing.
       between1 = round(tsat(1,n),4);
       between2= round(tsat(1,n+1),4);
       
       between0=between1-step;
       
       if round(fsat(1,n-1),4) == -1
          between0 = between2 + step;
       end
       if round(fsat(1,n-1),4)==1
          between0=between1-step;
       end
       
       if (round(fsat(1,n-1),4)==1 && round(fsat(1,n+1),4) ==1) || (round(fsat(1,n-1),4)==-1 && round(fsat(1,n+1),4) ==-1)
          between0=between1-step; 
       end
       % for ...... regression to get said points.
       for i=1: missing_points
              %linear regression using points between1 and between2.
            if round(between1 + (i*step),2) == 1.21 || (round(between1 + (i*step),2) >=2.09 && round(between1 + (i*step),2) <= 2.12)
             matA=[1, between1; 1, between2];
             matB=[round(fsat(1,n),4);round(fsat(1,n+1),4)];
             linear_values = matA\matB;
             
             regresion_y=(linear_values(1,1) + (linear_values(2,1)* (between1 + (i*step))));
             tsat_toadd = [tsat_toadd (between1 + (i*step))];
             fsat_toadd = [fsat_toadd (regresion_y)];
            else
                matA=[1, between0 + (i-1)*step , (between1^2); 1, between1 + ((i-1)*step),((between1 + ((i-1)*step))^2) ; 1, between2, (between2^2)];
              if between0>between2
                 matB=[round(fsat(1,n+2),4);round(fsat(1,n));round(fsat(1,n+1),4)];
              else
                 matB=[round(fsat(1,n-1),4);round(fsat(1,n));round(fsat(1,n+1),4)];
              end
              linear_values = matA\matB;
              regresion_y= linear_values(1,1) + (linear_values(2,1)* (between1 + (i*step))) + (linear_values(3,1)*(between1+(i*step))^2);
              tsat_toadd = [tsat_toadd (between1 + (i*step))];
              fsat_toadd = [fsat_toadd (regresion_y)];
                
            end  
       end
    end
end

%% Including calculated missing points in the original array
tsat = sort([tsat tsat_toadd]);
tsat(1)=[];
tsat=round(tsat, 4);
tsat_toadd=round(tsat_toadd, 4);
veces=0;
for x=2:1: size(tsat_toadd,2)
   
    for y=2: size(tsat,2)
       if tsat(1,y) ==  tsat_toadd(1,x)
           veces=veces+1;
          fsat= [ fsat(1:(y-1)) fsat_toadd(1,x) fsat(y:end)]; 
          break;
       end
    end
end

%% Local regression

%the calculated points are optimized to fit the original graph

pre=false;

previous_values=zeros(4,7); % reg 1 first x, reg 2 first f(x) reg 3 second x, reg 4 = f(second x)
post_values=zeros(4,7);

contador=1;
% we will need x^2, x^3 and x^4.

for x=1:180
    %we will scan for our "previous" and "post" values for every
    %crest/trough
    if (round(fsat(1,x),2)==1 || round(fsat(1,x),2)== -1) && pre==false
       %I'm already on the crest/trough. I need the previous 2 values. 
       previous_values(1,contador)=tsat(1,x-1);
       previous_values(2,contador)=fsat(1,x-1);
       previous_values(3,contador)=tsat(1,x-2);
       previous_values(4,contador)=fsat(1,x-2);
       pre=true;
    end
    if(round(fsat(1,x),2)<1 && round(fsat(1,x),2) > -1) && pre==true 
        %I can now store my post values.
        post_values(1,contador)=tsat(1,x);
        post_values(2,contador)=fsat(1,x);
        post_values(3,contador)=tsat(1,x+1);
        post_values(4,contador)=fsat(1,x+1);
        contador=contador+1;
        pre=false;
    end 
    
end
% we now have our regression points. Time to make the regression

%% Doing the general regression

% 4th Degree regression was used to obtain the equation of each parabola

n =4; % number of data points
mat_reg=zeros(5,5);
mat_sol=zeros(5,1);
results=0;
sumX=zeros(8,7);
sumY=zeros(5,7);
iterator_reg=1;
fun=zeros(7,5);
for x=1:length(post_values)
    for i=1:8
        if i<=5
            sumY(i,x)=sumatoryY(previous_values, post_values, x, i-1);
        end
     sumX(i,x)= sumatoryX(previous_values, post_values,x,i);
    end 
    clear mat_reg;
    clear mat_sol;
    
    %Solving the matrix
    mat_reg = [n, sumX(1,x), sumX(2,x), sumX(3,x), sumX(4,x);
               sumX(1,x), sumX(2,x), sumX(3,x), sumX(4,x), sumX(5,x);
               sumX(2,x), sumX(3,x), sumX(4,x), sumX(5,x), sumX(6,x);
               sumX(3,x), sumX(4,x), sumX(5,x), sumX(6,x), sumX(7,x);
               sumX(4,x), sumX(5,x), sumX(6,x), sumX(7,x),sumX(8,x)];
           
    mat_sol = [sumY(1,x);sumY(2,x);sumY(3,x);sumY(4,x);sumY(5,x)];
    mat_reg=vpa(mat_reg);
    mat_sol=vpa(mat_sol);
    solut= vpa((mat_reg\(mat_sol)));
    for t=1:5
       fun(x,t)=solut(t,1); 
    end
    iterator_reg=iterator_reg+1;
end

already=false;

%Improving precision
for x=1: 7
    for num=1:301
        if round(fsat(1,num),3)==1 || round (fsat(1,num),3) == -1
             fsat(1,num)= vpa(fun(x,1) + (fun(x,2)*tsat(1,num)) + (fun(x,3)*tsat(1,num)^2) + (fun(x,4)*tsat(1,num)^3) + (fun(x,5)*tsat(1,num)^4));
             already=true;
             if( round(fsat(1,num+1),3)<1 && round(fsat(1,num+1),3)> -1) && already==true
                 already=false;
                 break;
             end
        end
    end
end

plot(tsat,fsat,'*-');

%% Finding crests and valleys

%move along x in the range to find crests and valleys
crestsX = zeros(1,1);
crestsY= zeros(1,1);
valleyX = zeros(1,1);
valleyY= zeros(1,1);

%A critical point (maximum or minimum) is obtained if the previous and the
%next value from the iteration ara lower/higher than the current value
for i=2:180
   if (fsat(i) > fsat(i-1) && fsat(i)> fsat(i+1))
       if crestsX(1,1)==0
        crestsX(1,1)=tsat(i);
        crestsY(1,1)=fsat(i);
       else
        crestsX= [crestsX tsat(i)];
        crestsY=[crestsY fsat(i)];
       end
   end
   
   if (fsat(i) < fsat(i-1) && fsat(i)< fsat(i+1))
       if valleyX(1,1)==0
        valleyX(1,1)=tsat(i);
        valleyY(1,1)=fsat(i);
       else
        valleyX= [valleyX tsat(i)];
        valleyY=[valleyY fsat(i)];
       end
   end
end

disp("Crests: ");
for x=1:length(crestsX)
   disp(['Crest in x = ' num2str(crestsX(x))  ' with f(x) = ' num2str(crestsY(x))]);
end
disp(" ");
disp("Valleys: ");
for y=1: length(valleyY)
       disp(['Valley in x = ' num2str(valleyX(y))  ' with f(x) = ' num2str(valleyY(y))]);

end

%% Finding inflection points

diffFsat = backwardDiff(fsat, step);
inflectionpts = infPoints(diffFsat, tsat, fsat);
disp(" ");
disp("Inflection points");
for i=1: 5
    disp("x: " + inflectionpts(1,i) + " f(x): " + inflectionpts(2,i));
end

%% Finding times beyond the limits
times = countTime (fsat, step);
disp(" ");
for i=1:2
    if (i==1)
        disp("Time in which function is above 1.5: " + times(1,i));
    end
    if (i==2)
        disp("Time in which function is below -1: " + times(1,i));
    end
end


%% functions

%Gets the summatory in X for the matrix
function [a0]= sumatoryX(previous_values, post_values, col, pow)
a0=0;
for i=1:2:4
   a0 = a0 + (previous_values(i,col)^pow);
   a0= a0 + (post_values(i,col)^pow);
end
end

%Get the summatory in Y for the matrix
function [a0]= sumatoryY(previous_values, post_values, col, pow)
a0=0;
for i=2:2:4
   valuesx= previous_values(i-1,col);
   a0 = a0 + ((previous_values(i,col)) * (valuesx^pow) );
   valuesx= post_values(i-1,col);
   a0= a0 + (post_values(i,col) * (valuesx^pow));
end
end

%Returns the derivative of the function at a certain point
function [diffFsat] = backwardDiff(fsat, step)
diffFsat = zeros(1, 301);
for i=2:size(fsat,2)
    diffFsat(1, i) = (fsat(1,i)-fsat(1, i-1))/(step);
end
end

%Returns the inflection points
function [inflection2] = infPoints (diffFsat, tsat, fsat)
count = 1;
inflection = zeros(2, 301);
for i=2:size(diffFsat, 2)-1
    if(diffFsat(1, i-1) > diffFsat(1, i) && diffFsat(1,i+1) > diffFsat(1, i))
        inflection(1,count) = tsat(1, i);
        inflection(2, count) = fsat(1, i);
        count = count+1;
    end
    
    if(diffFsat(1, i-1) < diffFsat(1, i) && diffFsat(1, i+1) < diffFsat(1, i))
        inflection(1,count) = tsat(1, i);
        inflection(2, count) = fsat(1, i);
        count = count+1;
    end
end

count = 1;
inflection2 = zeros(2, 301);
for i=1:size(inflection, 2)-1
    if(inflection(1, i+1) - inflection(1, i) > 0.05)
        inflection2(1, count) = inflection(1, i);
        inflection2(2, count) = inflection(2, i);
        count = count + 1;
    end
end
end

%Returns the times beyond the limits
function [times] = countTime (fsat, step)
times = zeros(1,2);
time1 = 0;
time2 = 0;
for i=1:size(fsat, 2)
    if (fsat(1, i) > 1.5)
        time1 = time1 + step;
    end
    if (fsat(1, i) < -1)
        time2 = time2 + step;
    end
end
times(1,1) = time1;
times(1,2) = time2;
end