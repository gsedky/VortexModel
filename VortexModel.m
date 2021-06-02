function VortexModel(Opt, Duration, n, alpha_input, SurgeProfile, FinalSpeed, ...
                     Acceleration, PitchProfile, Pivot, omega, GustOpt, ...
                     GustWidth, GustType, C, X_initial, xmin, xmax, ymin,...
                     ymax, clmin, clmax, GR, Pivot_sin, omega_sin, ...
                     PitchAmplitude, Pivot_Control, K_p, RS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will serve as the main script to my vortex model code. The
% aim here is to be able to model the lift on a wing encountering a gust.

% Notes:
% The values with the subscript 'b' are in the body frame, values not
% with a subscript 'b' are in the inertial frame

% Inputs:
% Opt: TE shedding or LE and TE shedding
% Duration: Total time run for simulation
% n: number of control points on the wing
% alpha: angle of attack
% Surge Profile: Choose between a number of different surge motion profiles
% FinalSpeed: Maximum speed of surge

% Girguis Sedky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intialize Code
% Define Parameters
deg2rad = pi/180;                        % Change degrees to radians
alpha = alpha_input*deg2rad;                   % AoA in radians
% memory allocation
C_l = 0;                                 % Allocate memore to the coefficicient of lift
I_y_old = 0;

% Identify the initial vortex points and the control points on the wing chord
% Body Frame!!!!
d = C/(2*n);                             % Length between a vortex and a control point
Ratio = 2;
r_c = d/Ratio;
% note that the first and last vortex elements here are not bound but
% nascent
n_v = n+1;

% Control points
% Control point Position
X_b = d:2*d:C-d;                                % The x positons of all the control points
c.X_b{1} = X_b.';                               % Turn row to column vector
c.Y_b{1} = zeros(n,1);                          % The y positons of all the control points
% Transform control points to the inertial frame
[c.X{1}, c.Y{1}] = Trans_bi(c.X_b{1}, c.Y_b{1}, alpha);

% Bound vortices
% Bound vortices positions and strenght
X_b = 0:2*d:C;                                  % The x positons of all the control points
vb.X_b{1} = X_b.';                              % Turn row to column vector
vb.Y_b{1} = zeros(n_v,1);                       % The y positons of all the control points
vb.G{1} =  zeros(n_v,1);                        % Strength of vortices

% Transform bound vortex positions to the inertial frame
[vb.X{1}, vb.Y{1}] = Trans_bi(vb.X_b{1}, vb.Y_b{1}, alpha);

% Shed vortices
vs.X{1} = [];                               % The x positons of all the control points
vs.Y{1} = [];                               % The y positons of all the control points
vs.G{1}= [];                                % Strength of vortices


%% Create time vector
U_final = FinalSpeed;                        % Final speed in m/s

% define delta t based on the surge final speed
k = U_final;
%dt=0.01;
dt = (2*d/(k));
Time = 0:dt:Duration;                        % Time vector

%% Define motion
% Step Response
if strcmp('Step Response', SurgeProfile)
    [U,V] = StepResponse(FinalSpeed,Time);
elseif strcmp('Trapezoidal', SurgeProfile)
    [U,V] = Trapezoid(U_final,Time,Acceleration);
elseif strcmp('Linear Ramp', SurgeProfile)
    [U,V] = LinearRamp(U_final,Time,Acceleration);
else
    error('Motion profile not yet supported.')
end


%% Pitching Motion

% No pitching motion
if strcmp('None', PitchProfile)
    % Pitch axis
    a = Pivot*C;
    % inertial position vector of pitching point
    rp = [a*cos(alpha), -a*sin(alpha)];
    [c,vb,alpha] = NoPitching(alpha,Time,c,vb,U);
    
    % Pitching motion with a constant pitch rate
elseif strcmp('Constant', PitchProfile)
    % Pitch axis
    a = Pivot*C;
    % inertial position vector of pitching point
    rp = [a*cos(alpha), -a*sin(alpha)];
    [c,vb,alpha] = ConstantPitching(alpha,Time,c,vb,omega,a,Duration,U);
    
    % Pitching motion with sinusoid
elseif strcmp('Sinusoidal', PitchProfile)
    % Pitch axis
    a = Pivot_sin*C;
    % inertial position vector of pitching point
    rp = [a*cos(alpha), -a*sin(alpha)];
    [c,vb,alpha] = SinusoidalPitching(alpha,Time,c,vb, omega_sin,a,PitchAmplitude,U);

    % closed-loop pitch control scenario
elseif strcmp('Closed-loop', PitchProfile)
    % Pitch axis
    a = Pivot_Control*C;
    % inertial position vector of pitching point
    rp = [a*cos(alpha), -a*sin(alpha)];
%     [c,vb,alpha] = NoPitching(alpha,Time,c,vb,U);
    
    eps = 0;                 % Initialize eps
    C_l_ref = 0;             % regulate about 0 coefficient of lift
    omega_old = 0;           % initialize the pitch rate
    alpha_old = alpha*(pi/180);           % Initial AoA
end

%% Intialize circ
circ = 0;

% Create the matrix A that contains the coefficients of (gamma/2*pi*r)
% between bound vortices and control points, lamb oseen
Coeff = (-2*n+1:2:2*n-1);
row1= ((-1/(2*pi*d))./Coeff).*(1 - exp(-((d*Coeff).^2)/(r_c^2)));
A = gallery('circul',row1);
A = [A(1:n,n:end); ones(1,n_v)];

%% Intialize figures and animation
fig1 = figure(1);
% plot the wing with the vorticity and all that jazz
subplot(2,1,1);
h(1) = plot(vb.X{1},vb.Y{1},'k','LineWidth',2);
hold on;
h(2) = plot(vb.X{1}(1),vb.Y{1}(1),'.r','MarkerSize',10);
h(3) = plot(vb.X{1}(1),vb.Y{1}(1),'.b','MarkerSize',10);
xlabel('x, m','Fontsize',25,'interpreter','latex');
ylabel('y, m','Fontsize',25,'interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex','FontSize',20);
ylim([ymin ymax]);
xlim([xmin xmax]);
grid on;


% plot instantaeous lift evolution
subplot(2,1,2)
h(6) = plot(0,0,'.');
xlabel('Time, s','Fontsize',25,'interpreter','latex');
ylabel('$C_l$','Fontsize',25,'interpreter','latex');
xlim([0 Duration]);
ylim([clmin  clmax]);
set(gca,'TickLabelInterpreter', 'latex','FontSize',20);
grid on;
x0=100;
y0=-200;
width=3*550;
height=3*400;
set(fig1,'position',[x0,y0,width,height])
v = VideoWriter('Animation.avi');
open(v);

%% Start model
for ii=1:length(U)
    
        
    %% If closed-loop control is ON
    %% Control block
    % dont start input until leading edge hits the gust
    if strcmp('Closed-loop', PitchProfile)
        u(ii) = K_p*eps;
        %omega(ii) = 0.5*eps;
        omega(ii) = omega_old+u(ii)*dt;
        alpha(ii) = alpha_old+omega(ii)*dt;
        
        
        c.Vp_b{ii} = -omega(ii)*(c.X{1}-a);
        vb.Vp_b{ii} = -omega(ii)*(vb.X{1}-a);

    end
    
    

    
    %% Find the vertical velocity induced by the gust on control points , bound vortices and shed vorticity
    % This application of the gust velocity field is from AIAA 2019
    % journal paper
    
    % Initialize the vectors that will hold the relative velocity due to plunge
    % and gust at the control points, bound vortices and shed vortices
    V_cp = ones(size(c.X{ii}))*V(ii);    % Control points
    V_b = ones(size(vb.X{ii}))*V(ii);    % bound vortices
    V_s = ones(size(vs.X{ii}))*V(ii);    % Shed vortices
    % These are the terms with which you can figure out what parts of the wing
    % is in the boundary of the gust
    term = (U_final*Time(ii)-X_initial-c.X{ii})/GustWidth;      % Control points
    term_1 = (U_final*Time(ii)-X_initial-vb.X{ii})/GustWidth;   % bound vortices
    term_2 = (U_final*Time(ii)-X_initial-vs.X{ii})/GustWidth;   % Shed vortices
    if strcmp('On', GustOpt) % Gust option is turned on
        % Sine square gust
        if strcmp('Sine-square', GustType)
            % Control points
            V_cp(term>0 & term<1) = V_cp(term>0 & term<1) + GR*U_final*sin(term(term>0 & term<1)*pi).^2;
            % Bound vortices
            V_b(term_1>0 & term_1<1) = V_b(term_1>0 & term_1<1) + GR*U_final*sin(term_1(term_1>0 & term_1<1)*pi).^2;
            % shed vortices
            V_s(term_2>0 & term_2<1) = V_s(term_2>0 & term_2<1) + GR*U_final*sin(term_2(term_2>0 & term_2<1)*pi).^2;
            % Top hat gust
        elseif strcmp('Top-hat', GustType)
            % Control points
            V_cp(term>0 & term<1) = V_cp(term>0 & term<1) + GR*U_final;
            % bound vortices
            V_b(term_1>0 & term_1<1) = V_b(term_1>0 & term_1<1) + GR*U_final;
            % shed vortices
            V_s(term_2>0 & term_2<1) = V_s(term_2>0 & term_2<1) + GR*U_final;
        elseif strcmp('Trapezoid', GustType)
            % Control points
            index1=find(U_final*Time(ii)-X_initial-c.X{ii}<GR*U_final/RS & U_final*Time(ii)-X_initial-c.X{ii}>0);
            index2=find(U_final*Time(ii)-X_initial-c.X{ii}>GR*U_final/RS & U_final*Time(ii)-X_initial-c.X{ii}<GustWidth-GR*U_final/RS);
            index3=find(U_final*Time(ii)-X_initial-c.X{ii}>GustWidth-GR*U_final/RS & U_final*Time(ii)-X_initial-c.X{ii}<GustWidth);
            V_cp(index1) = V_cp(index1) + RS*(U_final*Time(ii)-X_initial-c.X{ii}(index1));
            V_cp(index2) = V_cp(index2) + GR*U_final;
            V_cp(index3) = V_cp(index3) + -RS*(U_final*Time(ii)-X_initial-c.X{ii}(index3)-GustWidth);
            
            % bound vortices
            index1=find(U_final*Time(ii)-X_initial-vb.X{ii}<GR*U_final/RS & U_final*Time(ii)-X_initial-vb.X{ii}>0);
            index2=find(U_final*Time(ii)-X_initial-vb.X{ii}>GR*U_final/RS & U_final*Time(ii)-X_initial-vb.X{ii}<GustWidth-GR*U_final/RS);
            index3=find(U_final*Time(ii)-X_initial-vb.X{ii}>GustWidth-GR*U_final/RS & U_final*Time(ii)-X_initial-vb.X{ii}<GustWidth);
            V_b(index1) = V_b(index1) + RS*(U_final*Time(ii)-X_initial-vb.X{ii}(index1));
            V_b(index2) = V_b(index2) + GR*U_final;
            V_b(index3) = V_b(index3) + -RS*(U_final*Time(ii)-X_initial-vb.X{ii}(index3)-GustWidth);
            % shed vortices
            index1=find(U_final*Time(ii)-X_initial-vs.X{ii}<GR*U_final/RS & U_final*Time(ii)-X_initial-vs.X{ii}>0);
            index2=find(U_final*Time(ii)-X_initial-vs.X{ii}>GR*U_final/RS & U_final*Time(ii)-X_initial-vs.X{ii}<GustWidth-GR*U_final/RS);
            index3=find(U_final*Time(ii)-X_initial-vs.X{ii}>GustWidth-GR*U_final/RS & U_final*Time(ii)-X_initial-vs.X{ii}<GustWidth);
            V_s(index1) = V_s(index1) + RS*(U_final*Time(ii)-X_initial-vs.X{ii}(index1));
            V_s(index2) = V_s(index2) + GR*U_final;
            V_s(index3) = V_s(index3) + -RS*(U_final*Time(ii)-X_initial-vs.X{ii}(index3)-GustWidth);
        end
    end
    %% Find the normal velocity at the control points in
    [~,c_v_b] = FindRelFlowVel(c.X{ii},c.Y{ii}, vs.X{ii}, vs.Y{ii}, vs.G{ii},r_c,U(ii),V_cp,alpha(ii),c.Vp_b{ii});
    
%     %% In this bit, I place nascent the vortex according to Ansari
%     if ii>1
%         [A,vb.X{ii},vb.Y{ii}]=AnsariPlacement(A,vb.X{ii},vb.Y{ii},c.X{ii},c.Y{ii},vs.X{ii},vs.Y{ii},r_c,alpha(ii),Opt);
%     end
    
    %% Find the new bound vortex sheet strength
    %obtain value for new bound vortex sheet strength
    vb.G{ii} = A\[-c_v_b; -circ];
    
    % Find the induced velocity at the free vortices due to all vortices
    % inertial frame, since it is in the inertial frame, you dont
    % incoprtate the velocity due to pitching
    [vs.u{ii},  vs.v{ii}] = velsLO(vs.X{ii},vs.Y{ii}, [vb.X{ii}; vs.X{ii}], [vb.Y{ii}; vs.Y{ii}], [vb.G{ii}; vs.G{ii}], r_c);
    % Find the induced velocity on the wing
    [vb.u{ii},  vb.v{ii}] = velsLO(vb.X{ii}, vb.Y{ii}, [vb.X{ii}; vs.X{ii}], [vb.Y{ii}; vs.Y{ii}], [vb.G{ii}; vs.G{ii}], r_c);
    
    
    % add the kinematic velocity constraint due to surge and plunge and gust to the shed and
    % bound vorticity
    % horizontal
    vs.u{ii} = vs.u{ii} + U(ii);
    vb.u{ii} = vb.u{ii} + U(ii);
    % vertical
    vs.v{ii} = vs.v{ii} + V_s;
    vb.v{ii} = vb.v{ii} + V_b;
    
    %% Update value for circulation in the flow
    circ = sum(vs.G{ii});
    %% Compute instantaneous lift
    
    % Compute impulse without density
    I_y = sum([vb.X{ii}; vs.X{ii}].*[vb.G{ii}; vs.G{ii}]);
    % Find the C_l using the impulse method
    L= (I_y-I_y_old)/dt;                 % Find the lift/rho by differentiating the impulse
    C_l = [C_l, 2*L./(C*U_final.^2)];
    % assign the current impulse to old impulse variable for the next loop
    I_y_old = I_y;
    
    
    
    %% animate solution and save recording
    % plot wing and vortices
    subplot(2,1,1);
    
    
    % Separate out the positive and the negative shed vortices
    % vortices with positive strength
    x_pos = vs.X{ii}(vs.G{ii}>0);
    y_pos = vs.Y{ii}(vs.G{ii}>0);
    % vortices with negative strength
    x_neg = vs.X{ii}(vs.G{ii}<0);
    y_neg = vs.Y{ii}(vs.G{ii}<0);
    
    if strcmp('On', GustOpt)
        % plot the gust patch
        eval('PlotPatch_FlowField');
    end
    
    set(h(1),'xData',c.X{ii},'yData',c.Y{ii});
    set(h(2),'xData',x_pos,'yData',y_pos);
    set(h(3),'xData',x_neg,'yData',y_neg);
    
 
    
    
    
    % plot the coefficient of lift
    subplot(2,1,2);
    set(h(6),'xData',Time(1:ii),'yData',C_l(1:ii));
    drawnow;
    M(ii) = getframe(fig1);
    writeVideo(v,M(ii));
    
    
%     %% Shed a vortex from edges and propagate shed vortices
%     % Enforce the kutta condition at the trailing and leading edge and convect wake
%     [vs.X{ii+1}, vs.Y{ii+1}, vs.G{ii+1}] = Shedding(vs.X{ii}, vs.Y{ii}, vs.G{ii}, vb.X{ii}, vb.Y{ii}, vb.G{ii}, vs.u{ii}, vs.v{ii}, vb.u{ii}, vb.v{ii}, dt, Opt);
    
    
    %% Shed a vortex from edges and propagate shed vortices
if strcmp('On', Opt)
    % Enforce the kutta condition at the trailing and leading edge and convect wake
    [vs.X{ii+1}, vs.Y{ii+1}, vs.G{ii+1}] = LE_Shedding(vs.X{ii}, vs.Y{ii}, vs.G{ii}, vb.X{ii}, vb.Y{ii}, vb.G{ii}, vs.u{ii}, vs.v{ii}, vb.u{ii}, vb.v{ii}, dt);
else
    % Enforce the kutta condition at the trailing edge and convect wake
    [vs.X{ii+1}, vs.Y{ii+1}, vs.G{ii+1}] = TE_Shedding(vs.X{ii}, vs.Y{ii}, vs.G{ii}, vb.X{ii}, vb.Y{ii}, vb.G{ii}, vs.u{ii}, vs.v{ii}, vb.u{ii}, vb.v{ii}, dt);
end


    %% Move the wing in the inertial frame (pitching)
    % Control points
    [rc_x, rc_y] = Trans_bi(c.X_b{ii}-a,0, alpha(ii));
    c.X{ii+1} = rp(1) + rc_x;
    c.Y{ii+1} = rp(2) + rc_y;
    c.X_b{ii+1} = c.X_b{ii};
    c.Y_b{ii+1} = c.Y_b{ii};
    % bound vortices
    [rc_x, rc_y] = Trans_bi(vb.X_b{ii}-a, 0, alpha(ii));
    vb.X{ii+1} = rp(1) + rc_x;
    vb.Y{ii+1} = rp(2) + rc_y;
    vb.X_b{ii+1} = vb.X_b{ii};
    vb.Y_b{ii+1} = vb.Y_b{ii};
    
    
    
    % if closed-loop control is on, update the lit error
    if strcmp('Closed-loop', PitchProfile)
    eps = C_l_ref-C_l(ii);
    omega_old = omega(ii);
    alpha_old = alpha(ii);
    end
end
close(v);

%% Save 
% Write options used into a text file
fileID = fopen('Parameters.txt','w');
TXT = 'Opt = %s \nDuration = %.2f \nn = %.2f \nalpha = %.2f \nSurgeProfile = %s \nFinalSpeed = %.2f \nAcceleration = %.2f \nPitchProfile = %s \nGustOpt = %s \nGustWidth = %.2f \nGustType = %s \nGR = %.2f \nC = %.2f \nPivot = %.2f \nomega = %.2f \nPivot_sin = %.2f \nomega_sin = %.2f \nPitchAmplitude = %.2f \nPivot_Control = %.2f \nK_p = %.2f';
fprintf(fileID,TXT,Opt,Duration,n,alpha_input,SurgeProfile,FinalSpeed,Acceleration,PitchProfile,GustOpt,GustWidth, GustType,GR,C,Pivot,omega,Pivot_sin,omega_sin,PitchAmplitude,Pivot_Control,K_p);


fclose(fileID);
% save the data
save('Data.mat','C_l','alpha','Time','C','U','V','vb','vs','c');
if strcmp('Closed-loop', PitchProfile)
    save('Data.mat','u', '-append')
end
close all;
end
