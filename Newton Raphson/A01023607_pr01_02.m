% Alberto Pascal A01023607
% Saúl E. Labra A01020725
% Rodrigo García A01024595
% Manuel Guadarrama A01020829
clc;
clear all;
close all;
Xi=1;
f= @(x) x^3 - (31/10)*x^2 + (1/10)*x + (21/5);
df = @(x) 3*x^2 - (31/5)*x + (1/10);


Xi_values=zeros(2,3);
Xi_values(1,1)=0.5;
Xi_values(1,2)=0.0161;
Xi_values(1,3)=2.051;
Xi_values(2,1)=0.5;
Xi_values(2,2)=0.0161;
Xi_values(2,3)=2.051;
solutions=zeros(2,3);
Err_tol=.000001;
initial_error=10000;
error = 1000000;
constant=0.5;
Err_vals=zeros(2,3);
Iterations=zeros(2,3);

%we calculate the normal values the newton raphson would give
for x=1:1:3
    Xi= Xi_values(1,x);
    disp(['First root for ' num2str(Xi) ' : ']);
    for n=0:1:100
        num1=f(Xi);
        num2=df(Xi);
        xnew= Xi- (num1/num2)*constant;
        if(n>=1)

            error= abs(((xnew-Xi)/xnew));
        end
        if (error<=Err_tol)
            Xi=xnew;
           break;
        end
        Xi=xnew;
       % disp ([' Mi Xi ahora va a valer ' num2str(Xi)]);

    end
    solutions(1,x)=Xi;
    Iterations(1,x)=n;
    Err_vals(1,x)=error;
    disp(['My approximation for X= ' num2str(Xi) ' with error of ' num2str(error) '% after iteration : ' num2str(n)]);
    error=10000;
end

for x=1:1:3
    Xi= Xi_values(1,x);
    disp(['Second root for ' num2str(Xi) ' : ']);
    for n=0:1:50000 % this number of cycles are because of constant2. We divided by 10,000 our original value. Therefore, 10,000 times more cycles.
        num1=f(Xi);
        num2=df(Xi);
        
        if n==0
           if num2>=0
               sign=-1;
           else
               sign=1;
           end
        end
        if(num2>0&& sign==-1) || (num2<0 && sign==1)
            constant=-0.5;
            constant2=-0.000154; % this value is obtained by doing Xmax-X where xmax= 0.016254 and x=0.0161
        else
            constant=0.5;
            constant2=0.000154;
        end
        if x==2
            xnew= Xi-(num1/num2)*constant2;
        else
            xnew= Xi - (num1/num2)*constant;
        end
        if(n>=1)
            
            error= abs(((xnew-Xi)/xnew));
            if n==1
                initial_error=error;
            end
            if error>initial_error
                disp('No other roots for this point');
                break;
            end
        end
        if (error<=Err_tol)
            Xi=xnew;
           break;
        end
        Xi=xnew;
      
       % disp ([' Mi Xi ahora va a valer ' num2str(Xi)]);

    end
    solutions(2,x)=Xi;
    Iterations(2,x)=n;
    Err_vals(2,x)=error;
    disp(['My approximation of X= ' num2str(Xi) ' with error of ' num2str(error) '% after iteration : ' num2str(n)]);
    error=10000;
end

% second part x=-2 and x=3.
Xi_2_values=zeros(2,2);
Xi_2_values(1,1)=-2;
Xi_2_values(2,1)=-2;
Xi_2_values(1,2)=3;
Xi_2_values(2,2)=3;
Err_vals_2=zeros(2,3);
Iterations_2=zeros(2,3);
solutions_2=zeros(2,3);
disp('Last point');


for x=1:1:2
    Xi= Xi_2_values(1,x);
    disp(['First root for ' num2str(Xi) ' : ']);
    for n=0:1:100
        num1=f(Xi);
        num2=df(Xi);
        xnew= Xi- (num1/num2)*constant;
        if(n>=1)

            error= abs(((xnew-Xi)/xnew));
        end
        if (error<=Err_tol)
            Xi=xnew;
           break;
        end
        Xi=xnew;
       % disp ([' Mi Xi ahora va a valer ' num2str(Xi)]);

    end
    solutions_2(1,x)=Xi;
    Iterations_2(1,x)=n;
    Err_vals_2(1,x)=error;
    disp(['My approximation for X= ' num2str(Xi) ' with error of ' num2str(error) '% after iteration : ' num2str(n)]);
    error=10000;
end

for x=1:1:2
    Xi= Xi_2_values(1,x);
    disp(['Second root for ' num2str(Xi) ' : ']);
    dont_print=false;
    for n=0:1:50000 % this number of cycles are because of constant2. We divided by 10,000 our original value. Therefore, 10,000 times more cycles.
        num1=f(Xi);
        num2=df(Xi);
        
        if n==0
           if num2>=0
               sign=-1;
           else
               sign=1;
           end
        end
        if(num2>0&& sign==-1) || (num2<0 && sign==1)
            constant=-0.5;
            constant2=-0.000154; % this value is obtained by doing Xmax-X where xmax= 0.016254 and x=0.0161
        else
            constant=0.5;
            constant2=0.000154;
        end
        if x==2
            xnew= Xi-(num1/num2)*constant2;
        else
            xnew= Xi - (num1/num2)*constant;
        end
        if(n>=1)
            
            error= abs(((xnew-Xi)/xnew));
            if n==1
                initial_error=error;
            end
            if error>initial_error
                disp('No other roots for this point');
                dont_print=true;
                break;
            end
        end
        if (error<=Err_tol)
            Xi=xnew;
           break;
        end
        Xi=xnew;
      
       % disp ([' Mi Xi ahora va a valer ' num2str(Xi)]);

    end
    solutions_2(2,x)=Xi;
    Iterations_2(2,x)=n;
    Err_vals_2(2,x)=error;
    if dont_print==false
        disp(['My approximation for X= ' num2str(Xi) ' with error of ' num2str(error) '% after iteration : ' num2str(n)]);
    end
    error=10000;
end

