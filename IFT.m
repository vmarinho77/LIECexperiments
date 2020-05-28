clear, close all, clc;

fname = 'IFT_experiment.json';
infos = jsondecode(fileread(fname));

%% Loading File
%General Informations
Ts = infos.Sample_time;
samples = infos.Samples;
iterations = infos.Iterations;

s = tf('s');
z = tf('z', Ts);

% Real Process   
numG = infos.numG;
denG = infos.denG';
G = tf(numG, denG, 'ioDelay', infos.Tau_dG); %G(s)
G = c2d(G, Ts);                              %G(z)

%Initial Controller
Kp_initial = infos.Kp_initial;
Ki_initial = infos.Ki_initial;
theta_initial = [Kp_initial Ki_initial];

%Reference Model
Am = infos.Am;           %Desired gain
beta = 2*Am/pi - 1;      %Beta
phi_m = pi/2*(1 - 1/Am); %Phase Margin
tau_d = infos.Tau_d;
tau_c = beta*tau_d;
Td = tf(1, [tau_c 1], 'ioDelay', tau_d);
Td = c2d(Td, Ts);

%% Initialization of IFT algortithm
time = 0: Ts: samples*Ts;   %Time simulation
time(length(time)) = [];

% time_period = 0: Ts: samples*Ts/2;
% time_period(length(time_period)) = [];

%Reference signal
r(:,1) = infos.Ref_amplitude/2 + infos.Ref_amplitude/2*square(pi*time/(infos.Ref_period*Ts));

yd = lsim(Td, r);           %Desired Output

%Class of Controller
PI_type = [1; Ts/(z-1)]; %Class of Controller = PI_Type
dC_dp = PI_type;         %C = p*PI_type => dC_dp = PI_type

%Iterative method
iterative_method = infos.Iterative_method;
init_step = infos.Initial_condition;    % 1 -> Steepest Descent 
                                        % 2 -> Aproximate Newton Raphson
%White Noise
noise_var = infos.Noise_var;

%Matrices initialization
dJ_dp = zeros(iterations,length(theta_initial));
dy_dp = zeros(samples,1);
J = zeros(1, iterations);
step_size = zeros(iterations,1);
theta = zeros(iterations + 1, length(theta_initial));
theta(1,:) = theta_initial;         %Initial controller parameters 

% figure
% y = lsim(feedback(p_initial*PI_type*G,1), r(:,1));
% yd = lsim(Td, r(:,1));
% plot(time, y), hold on
% plot(time, yd);


%% IFT algorithm
for i=1:iterations
    C(i) = theta(i,:)*PI_type;      %Transfer functuion od the Controller at each iteration
    T = feedback(C(i)*G, 1);    %Closed loop
    S = feedback(1, C(i)*G);    %Closed loop for the noise
    
    noise = sqrt(noise_var)*randn(1, samples); %Noise signal
    
    y(:,1) = lsim(T, r(:,1)) + lsim(S,noise);  %Output for each iteration
    u(:,1) = lsim(C(i), r(:,1) - y(:,1));      %Input for each iteration
    
    r(:,2) = r(:,1) - y(:,1);                  %Reference for the second experiment 
    y(:,2) = lsim(T,r(:,2)) + lsim(T,noise);   %Output of the second experiment 
    
    Q = 1/C(i)*dC_dp;                          %Auxiliar 
    
    dy_dp = lsim(Q,y(:,2));                     

    dJ_dp(i,:) = 2/samples *(y(:,1) - yd)'*dy_dp;
    
    J(i) = 1/samples*(y(:,1)-yd)'*(y(:,1)-yd);
    
    if iterative_method == 1
        if i>1
            if J(i)< J(i-1)
                step_factor = step_factor + 1;
            end
        end
        step_size(i) = init_step/step_factor;
    else
        R  = 1/samples*(dy_dp'*dy_dp); 
        step_size(i) = infos.Initial_condition;
    end
    theta(i+1, :) = theta(i, :) - step_size(i)*dJ_dp(i,:)*inv(R);  %Calculation of controller parameters for each iteration
end

%% Final Informations

%Saving Data
save('theta');
save('J');

C(i) = theta(end,:)*PI_type;  %Calculus of final C(z) 
T = feedback(C(i)*G, 1);      %Calculus of final T(z) 
S = feedback(1, C(i)*G);      %Calculus of final S(z) 

noise = sqrt(noise_var)*randn(1,samples);   %Noise signal
y(:,1) = lsim(T, r(:,1)) + lsim(S, noise);  %Final output 

final_cont = theta(end, :)   %Controller parameters after tuned
J_final = J(end)             %Cost Function after tuned

figure
hold on
y_initial(:,1) = lsim(feedback(G*C(1),1), r(:,1)) + lsim(feedback(1,G*C(1)),noise);
plot(time, y_initial(:,1));
plot(time,yd,'red');
plot(time,y(:,1));
legend('Saída inicial', 'Saída desejada','Saída adquirida',  'Location', 'southwest');
title('Resposta à referência utilizada no IFT');
xlabel('Tempo (segundos)'), ylabel('Amplitude');
grid
saveas(gcf, 'Saidas.png')

figure
plot(theta)
title('Evolução dos parâmetros do controlador a cada iteração IFT')
xlabel('Iterações')
legend('kp','ki','Location','northeast');
grid
saveas(gcf, 'Evolucao_theta.png')

figure
plot(J)
title('Evolução da função de custo a cada iteração IFT')
xlabel('Iterações')
saveas(gcf, 'Evolucao_J.png')
