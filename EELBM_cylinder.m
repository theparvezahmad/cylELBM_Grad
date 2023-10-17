% Entropic LBM for simulating flow past a circular cylinder. Free-slip
% boundary condition applied at the top and bottom walls. 

% The BGK collision operator is used. For entropic, refer to other two
% codes.

clc; clearvars; close all;

%% Initialize the variables
Re_d = 1000; tau = 0.51; u_inf = 0.085;
dia = 25; delta = 0.5; t_steps = 70000; viscosity = u_inf*dia/Re_d;
mulx = 25; muly = 18;
L = mulx*dia; H = muly*dia; nx = L+1; ny = H+1+2*delta; rad = 0.5*dia;
rho0 = 1.0; beta = 1.0/(1.0 + 6.0*viscosity);

f = zeros(ny,nx,9); ft = zeros(ny,nx,9); feq = zeros(9,1); wt=zeros(9,1); isn=zeros(ny,nx); kb=zeros(9,1);
vel_x = zeros(ny,nx); vel_y = zeros(ny,nx); vel_old = zeros(ny,nx); density = ones(ny,nx);
alpha = 2.0*ones(ny, nx); Cd = zeros(t_steps,1); Cl = zeros(t_steps,1);
Pxx = zeros(ny,nx); Pyy = zeros(ny,nx); Pxy = zeros(ny,nx);

vel_x_old(1:ny,1:nx) = u_inf; vel_y_old(1:ny,1:nx) = u_inf;

ex=[0,1,0,-1,0,1,-1,-1,1];
ey=[0,0,1,0,-1,1,1,-1,-1];

fx0 = 0; fy0 = 0;

%% Defining the fully developed velocity profile
x(1,1) = 0; x(nx,1) = L; x(2:nx-1,1) = 1:L-1;
y(1,1) = 0; y(ny,1) = H; y(2:ny-1,1) = 0.5:H-0.5;

[X,Y] = meshgrid(x,y);

% Bounceback directions and weights
kb(1)=1; kb(2)=4; kb(3)=5; kb(4)=2; kb(5)=3; kb(6)=8; kb(7)=9; kb(8)=6; kb(9)=7;
wt(1)=4/9; wt(2:5)=1/9; wt(6:9)=1/36;

%% Inlet velocity profile (uniform velocity at inlet)
% u_in = 6*u_avg*(H*y - y.^2)/(H^2);

%% Initializing f based on inlet velocity
center_x = 7*dia; center_y = 0.5*muly*dia; cx = center_x+1; cy = center_y+1+delta;
cy1=round(cy-1.5*rad); cy2=round(cy+1.5*rad); cx3=round(cx-1.5*rad); cx4=round(cx+1.5*rad);

for i = 1:ny 
    for j = 1:nx

        % Identifying solid nodes
        if i>=cy1 && i<=cy2 && j>=cx3 && j<cx4
            pt=sqrt((cx-j)^2+(cy-i)^2);
            if pt<=rad; isn(i,j)=1; end
        end

        if i==1 || i==ny; isn(i,j) = 1; end
        
        ux = u_inf;
        for a=1:9
            uc_alpha = ux*ex(a);
            f(i,j,a) = wt(a)*rho0*(1 + 3*uc_alpha + 4.5*uc_alpha^2 - 1.5*ux*ux);
        end

    end
end

%% File handle
file_Cd = ['Re_', num2str(Re_d), '_Dia_',num2str(dia),'_data_Cd_Cl.dat'];
handle1 = fopen(file_Cd, 'w');
fprintf(handle1, "variables = time, rho_avg, Cd, Cl\n\n");

% file_qi = ['Dia_',num2str(dia),'_qi_distance_data.txt'];
% qi_handle = fopen(file_qi, 'w');
% fprintf(qi_handle, "cx = %.2f, cy = %.2f\n", cx, cy);

% file_new = ['Re_', num2str(Re_d), 'Dia_',num2str(dia),'_tho_target_Data.txt'];
% filenew_handle = fopen(file_new, 'w');
% fprintf(filenew_handle, "cx = %.2f, cy = %.2f\n", cx, cy);

%% Time loop starts here.
xi = zeros(9,1);

for tim = 1:t_steps
    
    
    
    % Collision - essesntially entropic
    for i = 2:ny-1
        for j = 1:nx
            rho = 0; ux = 0; uy = 0;
            
            for a = 1:9
                rho = rho + f(i,j,a);
                ux = ux + f(i,j,a)*ex(a);
                uy = uy + f(i,j,a)*ey(a);
            end
            ux = ux/rho; uy = uy/rho;
            vel_x_old(i,j) = ux; vel_y_old(i,j) = uy; density(i,j) = rho;
            

%            Collision - SRT
%             u_sq = ux^2 + uy^2;
%             for a = 1:9
%                 uc_alpha = ux*ex(a) + uy*ey(a);
%                 feq = wt(a)*rho*(1 + 3*uc_alpha + 4.5*uc_alpha^2 - 1.5*u_sq);
% 
%                 ft(i,j,a) = f(i,j,a) - (f(i,j,a) - feq)/tau;
%             end

            % Evaluating equilibrium distribution function
            for a=1:9
                feq_x = (2 - sqrt(1+3*ux*ux))*((2*ux + sqrt(1+3*ux*ux))/(1-ux))^ex(a);
                feq_y = (2 - sqrt(1+3*uy*uy))*((2*uy + sqrt(1+3*uy*uy))/(1-uy))^ey(a);

                feq(a) = rho*wt(a)*feq_x*feq_y;
            end
            
            % Calculating alpha (analytical solution)
            for a = 1:9; xi(a) = feq(a)/f(i,j,a) - 1; end
            xi_max = max(abs(xi));

            if xi_max < 0.001
                alpha(i,j) = 2.0;
            else
                a1 = 0; b1 = 0; c1 = 0;
                for a = 1:9
                    if xi(a) < 0
                        a1 = a1 + f(i,j,a)*0.5*xi(a)^3;
                    end

                    b1 = b1 + f(i,j,a)*0.5*xi(a)^2;
                    c1 = c1 + f(i,j,a)*((2*xi(a)^2)/(2+xi(a)));
                end

                alpha_lower = 2*c1/(b1 + sqrt(b1^2 - 4*a1*c1));

                a2 = 0; b_cof1 = 0; b_cof2 = 0; c_cof = 0;
                for a = 1:9
                    if xi(a) < 0
                        a2 = a2 + (f(i,j,a)*(xi(a)^3)/6);
                    else
                        temp1 = 2/(4 + alpha_lower*xi(a)) + 1/(4 + 2*alpha_lower*xi(a)) + 2/(4 + 3*alpha_lower*xi(a));

                        b_cof2 = b_cof2 + f(i,j,a)*((2/15)*alpha_lower*beta^2*xi(a)^3)*temp1;
                    end

                    b_cof1 = b_cof1 + f(i,j,a)*(xi(a)^2)/2;

                    temp2 = 60*xi(a)^2 + 60*xi(a)^3 + 11*xi(a)^4;
                    temp3 = 60 + 90*xi(a) + 36*xi(a)^2 + 3*xi(a)^3;
                    c_cof = c_cof + f(i,j,a)*(temp2/temp3);
                end
                a2 = a2*beta^2;
                b_cof = b_cof1 - b_cof2;

                h_root = 2*c_cof/(b_cof + sqrt(b_cof^2 - 4*a2*c_cof));

                a_cof = 0;
                for a = 1:9
                    if xi(a) < 0
                        temp4 = (xi(a)^3)/6 - (h_root*beta*xi(a)^4)/12 + (h_root^2*beta^2*xi(a)^5)/20 - (h_root^3*beta^3*xi(a)^6)/5;
                        a_cof = a_cof + f(i,j,a)*temp4;
                    end
                end
                a_cof = beta^2*a_cof;

                alpha_higher = 2*c_cof/(b_cof + sqrt(b_cof^2 - 4*a_cof*c_cof));

                alpha(i,j) = alpha_higher;

                % ensuring that the populations remain positive
                xi_min = min(xi);
                alpha_max = (-1)/(beta*xi_min);
                if alpha_higher > alpha_max
                    alpha(i,j) = 0.5*(1 + alpha_max);
                end

            end         % if ends - alpha calculation


            % Collision
            for a=1:9
                ft(i,j,a) = f(i,j,a) - alpha(i,j)*beta*(f(i,j,a) - feq(a));
            end     %entropic collision ends here.



        end     % j ends - collision
    end     % i ends - collision

    % Calculating pressure tensors
    for i = 1:ny
        for j = 1:nx
            xx = 0; xy = 0; yy = 0;
            for a = 1:9
                xx = xx + f(i,j,a)*ex(a)*ex(a);
                yy = yy + f(i,j,a)*ey(a)*ey(a);
                xy = xy + f(i,j,a)*ex(a)*ey(a);
            end
            Pxx(i,j) = xx; Pyy(i,j) = yy; Pxy(i,j) = xy;

        end
    end

    
    
    vel_x_old(ny,:) = u_inf; vel_x_old(1,:) = u_inf; % vel_x(:,1) = u_in; 
    vel_y_old(ny,:) = 0; vel_y_old(1,:) = 0;


    %% Streaming of the populations
    for i = 2:ny-1
        for j = 1:nx
            for a=1:9
                ia = i + ey(a); ja = j + ex(a);
                if ja>=1 && ja <= nx
                    f(ia,ja,a) = ft(i,j,a);
                end
            end
        end
    end

    % Velocity after streaming
    rho_avg=0; count_rho=0;
   
    for i = 2:ny-1
        for j = 1:nx
            rho = 0; ux = 0; uy = 0;
            
            for a = 1:9
                rho = rho + f(i,j,a);
                ux = ux + f(i,j,a)*ex(a);
                uy = uy + f(i,j,a)*ey(a);
            end
            ux = ux/rho; uy = uy/rho;
            vel_x(i,j) = ux; vel_y(i,j) = uy; %density(i,j) = rho;
            rho_avg = rho_avg + rho; count_rho = count_rho + 1;
        end
    end

    rho_avg = rho_avg/count_rho;
    rho_in_avg = mean(density(2:ny-1, 1));
    vel_x(ny,:) = u_inf; vel_x(1,:) = u_inf; vel_y(ny,:) = 0; vel_y(1,:) = 0;

    

    %% Applying the boundary conditions

    % Inlet boundary: Specify f as the equilibrium distribution function based on inlet velocity
    j = 1;
    for i = 2:ny-1
        ux = u_inf;

        for a=1:9
            uc_alpha = ux*ex(a);
            fx_eq = (1 + 3*uc_alpha + 4.5*uc_alpha*uc_alpha - 1.5*ux*ux);
            f(i,j,a) = density(i,j)*wt(a)*fx_eq;
        end
    end

    % Outlet: Grad's approximation (Europhysics 2006 by Chikatamarla, Ansumali and Karlin).
    j = nx;
    for i = 2:ny-1
        
        rho_ij = density(i,j);
        Pab = [Pxx(i,j), Pxy(i,j); Pxy(i,j), Pyy(i,j)];
        rho_cs2 = [rho_ij/3, 0; 0, rho_ij/3];
        cs2 = [1.0/3.0, 0; 0, 1.0/3.0];

        for a = [4,7,8]
            ci_alpha = [ex(a)*ex(a), ex(a)*ey(a); ex(a)*ey(a), ey(a)*ey(a)];

            mat1 = Pab - rho_cs2;
            mat2 = ci_alpha - cs2;

            mat_sum = 0;
            for row = 1:2
                for col = 1:2
                    mat_sum = mat_sum + mat1(row,col)*mat2(row,col);
                end
            end

            jc_alpha = rho_ij*(vel_x_old(i,j)*ex(a) + vel_y_old(i,j)*ey(a));

            f(i,j,a) = wt(a)*(rho_ij + 3*jc_alpha + 4.5*mat_sum);

        end
        

    end     % outlet boundary - i ends.

%     Cylinder surface
    for i = cy1:cy2
        for j = cx3:cx4
            uwx = 0; uwy = 0;

            check_node = 0;
            for a = 1:9
                ia = i-ey(a); ja = j - ex(a);
                if isn(ia,ja)==1
                    check_node = check_node + 1;
                    break;
                end
            end
            
            if isn(i,j) == 0 && check_node == 1

                % Terget density and target velocity
                rho1 = 0; rho2 = 0; utx = 0; uty = 0; ut_count = 0;
                diff_dx = 0; diff_dy = 0;
                for a = 1:9
                    ia = i + ey(a); ja = j + ex(a);
                    is = i - ey(a); js = j - ex(a);
    
                    if isn(is, js) == 1
                        rho1 = rho1 + f(i,j,kb(a));
    
                        ufx = vel_x(ia, ja); ufy = vel_y(ia, ja);

                        % Calculate qi here
                        i_fluid = i; j_fluid = j; i_solid = is; j_solid = js;
                        dist_diff = 10;

                        while dist_diff > 0.0001
                            i_mid = 0.5*(i_fluid + i_solid);
                            j_mid = 0.5*(j_fluid + j_solid);
                            dist_mid = sqrt((i_mid - cy)^2 + (j_mid - cx)^2);

                            if dist_mid < rad
                                i_solid = i_mid; j_solid = j_mid;
                            else
                                i_fluid = i_mid; j_fluid = j_mid;
                            end

                            dist_diff = abs(rad - dist_mid);

                        end     % while ends

                        dist_solid = sqrt((i - is)^2 + (j - js)^2);
                        dist_surf = sqrt((i - i_mid)^2 + (j - j_mid)^2);

                        qi = dist_surf/dist_solid;

                        if i==is && j~=js; diff_dx = diff_dx + 1; distance_dx = qi; end
                        if i~=is && j==js; diff_dy = diff_dy + 1; distance_dy = qi; end

%                         if tim==100
%                             fprintf(qi_handle, "i = %d\t j = %d\t a = %d\t is = %d\t js = %d\t qi = %.4f\t diff_dx = %d\t diff_dy = %d\t distdx = %.4f\t dist_dy = %.4f\n", i,j,a,is,js,qi,diff_dx, diff_dy, distance_dx, distance_dy);
%                         end
                        
                        utx = utx + (uwx + qi*ufx)/(1 + qi);
                        uty = uty + (uwy + qi*ufy)/(1 + qi);
                        ut_count = ut_count + 1;

                    else
                        rho2 = rho2 + f(i,j,a);
                    end     %if isn == 0 ends
                end         % for a=1:9 ends
    
                rho_t = rho1 + rho2;
                utx = utx/ut_count; uty = uty/ut_count;

%                 vel_x_old = vel_x; vel_y_old = vel_y;
                
                if i>=cy && j>=cx
                    quad = 1;
                    if diff_dx == 1
                        dudx = (vel_x_old(i,j+1) - uwx)/(1+distance_dx);
                        dvdx = (vel_y_old(i,j+1) - uwy)/(1+distance_dx);
                        optiondx = 1;
                    else
                        dudx = (vel_x_old(i,j+1) - vel_x_old(i,j-1))/2;
                        dvdx = (vel_y_old(i,j+1) - vel_y_old(i,j-1))/2;
                        optiondx = 2;
                    end

                    if diff_dy == 1
                        dudy = (vel_x_old(i+1,j) - uwx)/(1+distance_dy);
                        dvdy = (vel_y_old(i+1,j) - uwy)/(1+distance_dy);
                        optiondy = 1;
                    else
                        dudy = (vel_x_old(i+1,j) - vel_x_old(i-1,j))/2;
                        dvdy = (vel_y_old(i+1,j) - vel_y_old(i-1,j))/2;
                        optiondy = 2;
                    end

                elseif i>=cy && j<cx
                    quad = 2;
                    if diff_dx == 1
                        dudx = (uwx - vel_x_old(i,j-1))/(1+distance_dx);
                        dvdx = (uwy - vel_y_old(i,j-1))/(1+distance_dx);
                        optiondx = 1;
                    else
                        dudx = (vel_x_old(i,j+1) - vel_x_old(i,j-1))/2;
                        dvdx = (vel_y_old(i,j+1) - vel_y_old(i,j-1))/2;
                        optiondx = 2;
                    end

                    if diff_dy == 1
                        dudy = (vel_x_old(i+1,j) - uwx)/(1+distance_dy);
                        dvdy = (vel_y_old(i+1,j) - uwy)/(1+distance_dy);
                        optiondy = 1;
                    else
                        dudy = (vel_x_old(i+1,j) - vel_x_old(i-1,j))/2;
                        dvdy = (vel_y_old(i+1,j) - vel_y_old(i-1,j))/2;
                        optiondy = 2;
                    end

                elseif i<cy && j>=cx
                    quad = 4;
                    if diff_dx == 1
                        dudx = (vel_x_old(i,j+1) - uwx)/(1+distance_dx);
                        dvdx = (vel_y_old(i,j+1) - uwy)/(1+distance_dx);
                        optiondx = 1;
                    else
                        dudx = (vel_x_old(i,j+1) - vel_x_old(i,j-1))/2;
                        dvdx = (vel_y_old(i,j+1) - vel_y_old(i,j-1))/2;
                        optiondx = 2;
                    end

                    if diff_dy == 1
                        dudy = (uwx - vel_x_old(i-1,j))/(1+distance_dy);
                        dvdy = (uwy - vel_y_old(i-1,j))/(1+distance_dy);
                        optiondy = 1;
                    else
                        dudy = (vel_x_old(i+1,j) - vel_x_old(i-1,j))/2;
                        dvdy = (vel_y_old(i+1,j) - vel_y_old(i-1,j))/2;
                        optiondy = 2;
                    end

                else
                    quad = 3;
                    if diff_dx == 1
                        dudx = (uwx - vel_x_old(i,j-1))/(1+distance_dx);
                        dvdx = (uwy - vel_y_old(i,j-1))/(1+distance_dx);
                        optiondx = 1;
                    else
                        dudx = (vel_x_old(i,j+1) - vel_x_old(i,j-1))/2;
                        dvdx = (vel_y_old(i,j+1) - vel_y_old(i,j-1))/2;
                        optiondx = 2;
                    end

                    if diff_dy == 1
                        dudy = (uwx - vel_x_old(i-1,j))/(1+distance_dy);
                        dydy = (uwx - vel_y_old(i-1,j))/(1+distance_dy);
                        optiondy = 1;
                    else
                        dudy = (vel_x_old(i+1,j) - vel_x_old(i-1,j))/2;
                        dvdy = (vel_y_old(i+1,j) - vel_y_old(i-1,j))/2;
                        optiondy = 2;
                    end

                end

%                 if tim==100
%                     fprintf(filenew_handle, "i = %d\t j = %d\t quad = %d\t optiondx = %d\t optiondy = %d\t distance_dx = %.4f\t distance_dy = %.4f\t dudx = %.4f\t dvdx = %.4f\t dudy = %.4f\t dvdy = %.4f\n", i,j,quad,optiondx,optiondy,distance_dx,distance_dy,dudx,dvdx,dudy,dvdy);
%                 end

                u_ab = rho_t*[utx*utx, utx*uty; utx*uty, uty*uty];
                cs2 = [1.0/3.0, 0; 0, 1.0/3];

                vel_mat = (rho_t/(6*beta))*[2*dudx, (dvdx+dudy); (dudy+dvdx), 2*dvdy];

                % Update the unknown population
                for a = 1:9
                    ia = i - ey(a); ja = j - ex(a);

                    if isn(ia, ja) == 1
                        ci_alpha = [ex(a)*ex(a), ex(a)*ey(a); ey(a)*ex(a), ey(a)*ey(a)];

                        mat1 = u_ab - vel_mat;
                        mat2 = ci_alpha - cs2;
                        mat_sum  =0;

                        for row = 1:2
                            for col = 1:2
                                mat_sum = mat_sum + mat1(row,col)*mat2(row,col);
                            end
                        end

                        jc_alpha = rho_t*(ex(a)*utx + uty*ey(a));

                        f(i,j,a) = wt(a)*(rho_t + 3*jc_alpha + 4.5*mat_sum);
                    end

                end
            end

        end
    end


    i = 2;              % Bottom boundary - Free slip
    for j = 2:nx-1
        f(i,j,3) = f(i-1,j,5);
        f(i,j,6) = f(i-1,j,9);
        f(i,j,7) = f(i-1,j,8);
    end
    f(i,nx,3) = f(i-1,nx,5);

    i = ny-1;              % Top boundary - Free slip
    for j = 2:nx-1
        f(i,j,5) = f(i+1,j,3);
        f(i,j,8) = f(i+1,j,7);
        f(i,j,9) = f(i+1,j,6);
    end
    f(i,nx,5) = f(i+1,nx,3);
 

    %% Calculating forces (new)
   
    fx1 = 0; fy1= 0;
    for i = cy1:cy2
        for j = cx3:cx4
            if isn(i,j) == 0
                for a = 1:9
                    ia = i - ey(a); ja = j - ex(a);
                    if isn(ia, ja) == 1
                        fx1 = fx1 + ex(kb(a))*(f(i,j,a) + ft(i,j,kb(a)));
                        fy1 = fy1 + ey(kb(a))*(f(i,j,a) + ft(i,j,kb(a)));
                    end
                end
            end
        end
    end

    fx_avg = 0.5*(fx0 + fx1); fy_avg = 0.5*(fy0 + fy1);
    fy0 = fy1; fx0 = fx1;

    % Drag coefficients
    Cd(tim,1) = 2*fx_avg/(rho_in_avg*u_inf*u_inf*dia);
    Cl(tim,1) = 2*fy_avg/(rho_in_avg*u_inf*u_inf*dia);
    

    if mod(tim, 100) == 0
        fprintf(handle1, "%d\t%8.4f\t%8.4f\t%8.4f\n", tim, rho_avg, Cd(tim,1), Cl(tim,1));
    end

    if mod(tim,10000) == 0
        handle2 = fopen(['Re_', num2str(Re_d), '_Dia_',num2str(dia), '_Data_velocity_density', '.dat'],'w');
        fprintf(handle2,'variables=x,y,ux,uy,rho_avg\nZONE I=%d J=%d\n\n',nx,ny);
        for i = 1:ny
            for j = 1:nx
                fprintf(handle2,"%8.4f\t%8.4f\t%12.8f\t%12.8f\t%6.4f\n", x(j), y(i), vel_x_old(i,j), vel_y_old(i,j), density(i,j));
            end
        end
    end

%     vel_x_old = vel_x; vel_y_old = vel_y;
    fprintf("%d\n", tim);
end     % time loop ends here.

%% Writind the data file
filename = ['Re_', num2str(Re_d), '_Dia_',num2str(dia), '_Data_velocity_density', '.txt'];
handle2 = fopen(filename,'w');
fprintf(handle2,'variables=x,y,ux,uy,rho_avg\nZONE I=%d J=%d\n\n',nx,ny);
for i = 1:ny
    for j = 1:nx
        fprintf(handle2,"%8.4f\t%8.4f\t%12.8f\t%12.8f\t%6.4f\n", x(j), y(i), vel_x_old(i,j), vel_y_old(i,j), density(i,j));
    end
end

fclose('all');
contourf(X,Y,sqrt(vel_x.^2 + vel_y.^2));
set(gcf, 'position',[300,200,nx,ny]);
