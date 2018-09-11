% Alberto Pascal A01023607
% Saúl E. Labra A01020725
% Rodrigo García A01024595
% Manuel Guadarrama A01020829

clc;
clear all;
close all;

%Original function and its first and second derivative
f = @(x) x.^4 - 10.*x.^3 + 27.*x.^2 - 2.*x - 40;
df = @(x) 4.*x.^3 - 30.*x.^2 + 54.*x - 2;


%Error tolerance for NR method and initial error before execution
err_tol = 0.0001;
error = 100;

%Definition of the interval to find roots (-3,4)
lowest=-3;
highest=4;
%Stores the values where the function is going to be evaluated (intervals of 1)
values = lowest-1:1:highest;
%Stores the initial points for Newton Raphson
Xnew_arr=zeros(1,4);
%The first X for Newton Raphson will be the lower limit
Xnew_arr(1,1)=lowest;

%Solution stores the roots
solutions=zeros(1,4);
count_sol=1;
count=2;
first=true;
next=true;

%Cycle moves along the function in intervals of 1 from 2 (because code
%requires a "n-1" position in the array values) till the end of the array
%values
for n=2:1:length(values)
   
    %stores the previous and current iteration derivative of the function
    store_prev = df(values(n-1));
    store_curr=df(values(n));
    
    %In case there isn't a change of sign
    if (store_prev<=0 && store_curr<=0)
    end
    if (store_prev >=0 && store_curr>=0)
    end
    %In case there's a change of sign and it is the first value of the
    %interval
    if(store_prev >=0 && store_curr<=0) || (store_prev <=0 && store_curr>=0) && first==false
        first=true;
        next=false;
    end
    %In case there's a change of sign and is the first point of the
    %interval
    if(store_prev >=0 && store_curr<=0) || (store_prev <=0 && store_curr>=0&& next==false)
       
        %Calulates the point to evaluate Newton Raphson with an average of
        %the interval
        num_add=((values(n-1) + values(n))/2) ;
        
        %If the new value is lower than the lower limit of the interval
        if(num_add<lowest)
            num_add=lowest;
        end
        
        %If the new value is greater than the upper limit of the interval
        if(num_add>highest)
            num_add=highest;
        end
        
        %Stores the position in X where Newton Raphson is going to be
        %calculated
        Xnew_arr(count)= num_add;
        count=count+1;
    
    end
        
end


%Get the first root fo the function
%Cycle with 20 as the limit of iterations
for x=1:1:count-1
    Xi=Xnew_arr(x);
    for i=0:1:20
        num1 = f(Xi);
        num2 = df(Xi);
        Xnew = Xi - (num1/num2);

        %After the first iteration rel. Error can be calculated
        if(i>=1)
            error = abs((Xnew-Xi)/Xnew);
        end

        %When the error fulfills the tolerance, the cycle breaks
        if(error<=err_tol)
            break;
        end

        %The new value of Xi is assigned
        Xi = Xnew;
        disp (['It.#: ' num2str(i) ' Xi = ' num2str(Xi)]);
    end
    %The value of the root is dispalyed
    error=10000;
    solutions(count_sol)=Xnew;
    count_sol=count_sol+1;
    disp (['The root ' num2str(x) ' is ' num2str(Xnew)]);
end



