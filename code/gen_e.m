function [ x_arr,y_arr, vx_arr,vy_arr ] = gen_e( num, x_dim, y_dim, question )
%Generates the initial positions and velocities for questions 1,2,3

if question ==1
    %positions
    x_arr = rand(1,num);
    x_arr = x_arr*x_dim;
    y_arr = rand(1,num);
    y_arr = y_arr*y_dim;
    
    %velocity
    angle=rand(1,num)*2*pi;
    vx_arr = 132.2e3*cos(angle);
    vy_arr =132.2e3 * sin(angle);
    %vx_arr = (132.2e3)*rand(1,num);
    %vy_arr = sqrt(abs(132.2e3^2-vx_arr.^2));

    
elseif question == 2
    %velocity
    vx_arr = (132.2e3)*randn(1,num);
    vy_arr = (132.2e3)*randn(1,num);
    %position
    x_arr = rand(1,num);
    x_arr = x_arr*x_dim ;
    y_arr = rand(1,num);
    y_arr = y_arr*y_dim ;
    
    
elseif question ==3
    %velocities
    vx_arr = (132.2e3)*randn(1,num);
    vy_arr = (132.2e3)*randn(1,num);
    %positions
    x_arr = rand(1,num);
    x_arr = x_arr*x_dim ;
    y_arr = rand(1,num);
    y_arr = y_arr*y_dim ;

    %ensure no positions are inside the boundary
    q=1;
    while q<=length(x_arr)
        if y_arr(q)>0.6*y_dim && x_arr(q)>0.4*x_dim && x_arr(q)<0.6*x_dim
            y_arr(q) = rand(1)*y_dim;
            x_arr(q) = rand(1)*x_dim;
            q=q-1;
        end

        if q<1
            q=1;
        end

        if y_arr(q)<0.4*y_dim && x_arr(q)>0.4*x_dim && x_arr(q)<0.6*x_dim
            y_arr(q) = rand(1)*y_dim;
            x_arr(q) = rand(1)*x_dim;
            q=q-1;
        end       
        q=q+1;
    end
end
end 

