function Un = RK4(residual,t0,tf,h,Un,mesh)

% Description: Integrate ODE using RK4

az = 0; el = 0; % Used for plotting

t_next = t0 + h;
last = false;
while t_next <= tf
    fprintf('Next time: %f. Final time: %f.\n',t_next,tf)

    Rn = residual(Un,(t_next-h));
    
    U_hat = Un     + .5*h*Rn;
    R_hat = residual(U_hat,(t_next-h/2));
    
    U_wav = Un + .5*h*R_hat;
    
    R_wav = residual(U_wav,(t_next-h/2));
    
    U_bar     = Un + h*R_wav;
    R_bar = residual(U_bar,t_next);
    
    Un = Un + 1/6*h*(Rn + 2*(R_hat + R_wav) + R_bar); 
    
    %figure(2),plot_field(mesh,Un),  %view(2) %view([az,el]) 
    %drawnow
    %figure(3),plot_field(mesh,eqn.u(xxyy,t_next)), view([az,el])
    %drawnow

    % Get ready for next iteration
    if (t_next == tf)
        break;
    end
    t_next = t_next + h;
    if ~(last) && (t_next > tf) % Readjust time step if t_next > tf
        t_now = t_next - h;
        h = tf - t_now;
        t_next = t_now + h;

        last = true;
    end

end
    
end