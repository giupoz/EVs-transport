clear all
N_sim = 5;
plot_switch = 'ON';
save_data = 'OFF';

%% CONSTANTS
T_tot = 1200;              % total simulation time [s]
dt = 1e-4;                 % time step [s]
N_T = T_tot/dt;            % number of time steps
N = T_tot*2;               % size of x (I want the results every 0.5s)
k_B = 1.380649e-23;        % Boltzman's constant [J/K]
T = 310.15;                % temperature [K]
zeta = 0.0069;             % cytosol dynamic viscosity at 37°C [kg/(m*s)]
eta = 0.0006913;           % water dynamic viscosity at 37°C [kg/(m*s)]
R_r = 7e-9;                % receptor radius [m]
R = 300e-9;                % EV radius [m]
xi = 12*pi*R_r*zeta;       % friction coefficient of the receptor
gamma = 12*pi*R*eta;       % friction coefficient of the EV
xi_eff = (xi+gamma)/2;
k = 8e-9;                  % spring constant [N/m]

% V1 (asymmetric ratchet potential)
L1 = 5e-6;                 % period [m]
alpha1 = 1/5;              % symmetry factor
h1 = 4e+6*k_B*T;           % heigth
dV1_l = h1/alpha1;         % left derivative
dV1_r = -h1/(1-alpha1);    % right derivative
dV1 = 0;

% V2 (symmetric ratchet potential)
L2 = 2.3*R;                % period [m]
alpha2 = 1/2;              % symmetry factor% height
h2 = 5e+6*k_B*T;           % heigth
dV2_l = h2/alpha2;         % left derivative
dV2_r = -h2/(1-alpha2);    % right derivative
dV2 = 0;                   % current derivative

% transition rates
nu_on1 = 180;
nu_off1 = 500;

nu_on2 = 600;
nu_off2 = 10;

timev=[0:0.5:T_tot-0.5]';

for set=2:2 %CTRL=1, CYTOEV=2, CYTOHN=3
    

    if set==2
        beta = 1;
        nu_on = 0;
        nu_off = 0;
        folder = 'simulazioni/cytoEV_r3/';
    elseif set==3
        beta = 0;
        nu_on = 0;
        nu_off = 0;
        folder = 'simulazioni/cytoHN_r3/';
    else %set=1
        beta = 0;
        nu_on = 50;
        nu_off = 2400;
        folder = 'simulazioni/ctrl_r3/';
    end


 

    for n=1:N_sim 

        fprintf('set %d, simulation %d\n',set,n)

        pii=0;


        %unknowns
        x = zeros(N_T,1);
        xr = zeros(N_T,1);
        x_plot = zeros(N,1);
        xr_plot = zeros(N,1);
        states = beta*ones(N,1);
        j=1;

        for i=1:N_T-1

            if beta == 1 

                if pii == 0
                    
                    xr(i+1) = xr(i) + sqrt(k_B*T*dt/xi)*randn(1);
                    if rand < nu_on1*dt
                        pii=1;
                    end
                else %pii=1
                    if mod(xr(i),L1) >= alpha1*L1
                        dV1 = dV1_r;
                    else
                        dV1 = dV1_l;
                    end
                    
                    xr(i+1) = xr(i) - pii*dV1*dt/xi + sqrt(k_B*T*dt/xi)*randn(1);
                    if rand < nu_off1*dt
                        pii=0;
                    end
                end

                x(i+1) = x(i) - (k/gamma)*(x(i)-xr(i+1))*dt + sqrt(k_B*T*dt/gamma)*randn(1);

                if rand < nu_off*dt
                    beta = 0;
                end

            else %beta=0 ROLLING

                if pii == 0
                    
                    xr(i+1) = xr(i) + sqrt(k_B*T*dt/xi_eff)*randn(1);
                    if rand < nu_on2*dt
                        pii=1;
                    end
                else %pii=1
                    if mod(xr(i),L2) >= alpha2*L2
                        dV2 = dV2_r;
                    else
                        dV2 = dV2_l;
                    end
                   
                    xr(i+1) = xr(i) - pii*dV2*dt/xi_eff + sqrt(k_B*T*dt/xi_eff)*randn(1);
                    if rand < nu_off2*dt
                        pii=0;
                    end
                end

                x(i+1) = x(i) - (k/gamma)*(x(i)-xr(i+1))*dt + sqrt(k_B*T*dt/gamma)*randn(1);

                if rand < nu_on*dt
                    beta = 1;
                end

            end


            if mod(i*dt,0.5)==0

                j=j+1;
                states(j) = beta;
                x_plot(j) = x(i+1);
                xr_plot(j) = xr(i+1);
            end


        end

        x_plot_u=1e+6*x_plot;
        xr_plot_u=1e+6*xr_plot;
        if strcmp(plot_switch,'ON')

            figure(n)
            plot(timev,x_plot_u)

        end

        if strcmp(save_data,'ON')
            A=[x_plot_u, xr_plot_u, states];
            fileID = fopen([folder,sprintf('sim_%d',n),'.txt'],'w');
            fprintf(fileID, '%6s %12s %18s\n', 'x','xr','state');
            fprintf(fileID, '%+6f %+12f %14d\n',A');
            fclose(fileID);
        end

    end
end

