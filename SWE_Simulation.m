clear;clc;
t_end = 0.5;
% Grid creation
N = 60;
dx = 2/N;
dy = 2/N;
% x vector plus ghost cells
x = -1-dx:dx:1+dx;
y = x;
[xgrid,ygrid] = meshgrid(x,y);

% Initial Conditions
h = ones(size(xgrid));
for i = 1:length(x)
    for j = 1:length(y)
        if xgrid(i,j) >= -0.5 && xgrid(i,j) <= 0.5 && ygrid(i,j) >= -0.5 && ygrid(i,j) <= 0.5
             h(i,j) = 2;
        end
    end
end

%Unknown matrix M and matrices u,v
M_old = zeros([size(h) 3]); 
M_old(:,:,1) = h;
u_old = zeros(size(xgrid)); 
v_old = u_old;

h_int = [];

% Main iteration loop
time = 0;
while time < t_end
    
    % calculate lambda = abs(u or v) + sqrt(gh) => we use g = 1
    lambdax = abs(u_old(:,[2:63,1])) + sqrt(M_old(:,:,1));
    lambday = abs(v_old([2:63,1],:)) + sqrt(M_old(:,:,1));
    lammax = norm([lambdax(:); lambday(:)],Inf);
    
    dt = 0.5*0.5*(dx/lammax);
    time = time+dt;
    
    %Calculate f(v_i's) for lax freidrich
    % d/dx(hu,hu^2+gh^2/2,huv)
    comp1 = M_old(:,:,2);
    comp2 = M_old(:,:,2).^2./M_old(:,:,1)+0.5*M_old(:,:,1).^2;
    comp3 = M_old(:,:,2).*M_old(:,:,3)./M_old(:,:,1);
    ddx = cat(3,comp1,comp2,comp3);
    % d/dy(hv,huv,hv^2+gh^2/2)
    comp1 = M_old(:,:,3);
    comp2 = M_old(:,:,2).*M_old(:,:,3)./M_old(:,:,1);
    comp3 = M_old(:,:,3).^2./M_old(:,:,1)+0.5*M_old(:,:,1).^2;
    ddy = cat(3,comp1,comp2,comp3);

    x_flux =  0.5*(ddx+ddx(:,[2:63,1],:)) - 0.5.*lambdax.*(M_old(:,[2:63,1],:)-M_old);
    y_flux =  0.5*(ddy+ddy([2:63,1],:,:)) - 0.5.*lambday.*(M_old([2:63,1],:,:)-M_old);
    
    M_new = M_old - (dt/dx)*(x_flux - x_flux(:,[63,1:62],:)) - (dt/dy)*(y_flux - y_flux([63,1:62],:,:));
    
    % Boundaries - set edges by copying values
    % h
    M_old(63,1:63,1) =  M_old(62,1:63,1); 
    M_old(1,1:63,1) =  M_old(2,1:63,1);
    M_old(1:63,63,1) =  M_old(1:63,62,1); 
    M_old(1:63,1,1) =  M_old(1:63,2,1);
    % hu
    M_old(1:63,63,2) = -M_old(1:63,62,2); 
    M_old(1:63,1,2) = -M_old(1:63,2,2);
    M_old(63,1:63,2) =  M_old(62,1:63,2); 
    M_old(1,1:63,2) =  M_old(2,1:63,2);
    % hv
    M_old(1:63,63,3) =  M_old(1:63,62,3); 
    M_old(1:63,1,3) =  M_old(1:63,2,3);
    M_old(63,1:63,3) = -M_old(62,1:63,3); 
    M_old(1,1:63,3) = -M_old(2,1:63,3);

    % Plot
    surf(x,y,M_old(:,:,1))
    axis([-1 1 -1 1 0 3])
    xlabel x 
    ylabel y
    zlabel z
    pause(0.001)
    
    u_new = M_old(:,:,2)./M_old(:,:,1);   
    v_new = M_old(:,:,3)./M_old(:,:,1);
    M_old = M_new;
    u_old = u_new;
    v_old = v_new;
    
    % Variables for plot of h integrated over spatial domain
    total = 0;
    for i = 1:63
        for j = 1:63
            total = total + M_new(i,j,1)*dx*dy; 
        end
    end
    h_int = [h_int, total];
end


% Also for plot of h integrated over domain
%{
time_axis = [0];
for i = 2:length(h_int)
    time_axis = [time_axis, time_axis(i-1)+dt];
end

plot(time_axis,h_int)
axis([0 3 0 6])
title("h integrated over domain")
xlabel("Time")
ylabel("Integral of h")
%}
