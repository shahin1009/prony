clc
clear 
close all

%% Load H1 Data from each excitation
fsamp=10240; % Sampling frequency in Hz
df=0.1;
freq=0:df:fsamp/2;%nyquist

[f, p]=uigetfile('*.mat','Please specify the file to open');
cd(p);
load(f);

%% Plotting H1
figure

semilogy(freq,abs(H1))
axis tight
grid on
axis tight
xlabel('Frequency [Hz]')
ylabel('|H_1| [m/s^2/N]')
legend('N1', 'N2', 'N3','N4','N5')
grid on


%% MDOF Method
H1_10 = [H1(1,:);H1(2:end,:)/2;fliplr(conj(H1(2:end,:)/2))];

%% IRF

irf1 = ifft(H1_10, 'symmetric');

%% Main Prony method
dt=1/fsamp;
m=400;
order=100;
poles_num=8;

matrix1 = Prony(irf1,fsamp,order,m,H1,freq);
params_matrix1 =ParamFilter(matrix1,fsamp,poles_num);


%%
mode_number=5;
for i = 1:5
figure;
plot(ComplexModeToRealMode(params_matrix1(3:end,i)))
title('Mode shape at freq=',[num2str(params_matrix1(1,i))])
end
%% Normalize and real(ize) the mode shapes

ComplexModeToRealMode(params_matrix1(3:end,1));

%% the required hankel function to make the frequency histogram, also from another group

function [R,hhat] = hankel(IRF,N,m)
% Creating the Hmn(t) and Hmn(t+dt) matrix
ii = size(IRF,2);

G = [];
for k=1:ii
    irf = IRF(:,k);
    A1=zeros(m,1);
    A1(1:m,1) = irf(1:m);

    R1 = zeros(m,2*N);

    R1(:,1)=A1;

    for l = 1:2*N-1
        A1 = [A1(2:end,:);irf(m+l,1)];

        R1(:,l+1)=A1;
    end
    G = [G;R1];
end

R = G(:,1:2*N-1);
hhat = G(:,2*N);
end

function A = V_matrix(roots,m)
% Creating the Hmn(t) and Hmn(t+dt) matrix

A = zeros(m,length(roots));
for k=1:m
    A(k,:)=roots.^(k-1);
end

end

function matrix = Prony(IRF, fsamp, order, m, H1_1,freq)
    dt = 1 / fsamp;
    eigen_freqs = [];
    eigen_damps = [];
    sai=[];
    figure(10);
    % Plot H1_1 on the left y-axis
    yyaxis left
    semilogy(freq, abs(H1_1(:,1)))
    axis tight
    grid on
    xlabel('Frequency [Hz]')
    ylabel('|H_1| [m/s^2/N]')
    xlim([0,2500])
    % Initialize the second y-axis
    yyaxis right
    ylabel('Eigenfrequencies [Hz]')
    hold on
    
    % Prony method to compute eigenfrequencies and damping ratios
    disp('Doing it for the order of the sys:')
    for i = 2:order
        N = 2*i;
        if mod(i, 2) == 0
            disp(['N = ', num2str(i)]);
        end

        [R, hhat] = hankel(IRF, N, m);
        Beta = -pinv(R) * hhat;

        Beta = Beta(end:-1:1) ;

        Beta = [1; Beta] ;

        roots_v = roots(Beta);
        A = zeros(size(IRF,2),length(roots_v(:).'));


        for ii=1:size(IRF,2)

            A1=pinv(V_matrix(roots_v(:).',m))*R((ii-1)*m+1:ii*m,1);
            A(ii,:) = A1;
            
        end
        
        sai = [sai,A]; 

        pole_sr = log(roots_v) / dt;
        cisi = sqrt(1 ./ ((imag(pole_sr) ./ real(pole_sr)) .^ 2 + 1));
        eigen_damps = [eigen_damps, cisi'];  % Unrefined damping ratios
        
        freque = -real(pole_sr) ./ cisi / 2 / pi;
        eigen_freqs = [eigen_freqs, freque'];
        
        % Scatter plot of eigenfrequencies
        scatter(freque, ones(size(freque)) * i, 'filled')
    end

    hold off

    matrix = [eigen_freqs; eigen_damps;sai];
    
    % Plotting the frequency stability histogram
    freq_interval = 0.002 * fsamp / 2;  % based on 0.1% change
    freq_bins = 10:freq_interval:2500;

    figure;
    histogram(eigen_freqs, freq_bins, 'Edgecolor', 'k');
    xlabel('Frequency ranges');
    ylabel('Count');
    title('Frequency stability diagram');
end

function cond_matrix2 = ParamFilter(matrix,fsamp,pole_nums)


freq_interval = 0.002*fsamp/2;                    % based on 0.1% change
freq_bins = 10:freq_interval:2500;
[num1,edge1] = histcounts(matrix(1,:), freq_bins);
[~, indices] = sort(num1, 'descend');


edge1 = sort(edge1(indices(1:pole_nums)),'ascend');
cond_matrix=[];
cond_matrix2=[];
disp('The ranges of stable frequencies are:')

for i=1:pole_nums
    freq_Lower_bound = edge1(i);
    freq_Uper_bound = edge1(i)+freq_interval;


    disp([num2str(freq_Lower_bound), ...
        ' & ', ...
        num2str(freq_Uper_bound)] ...
        );


    % Specify the conditions for each row
    row1_condition = matrix(1, :) >= freq_Lower_bound & ...
        matrix(1, :) <= freq_Uper_bound;


    % Find the columns that satisfy all conditions
    columns_satisfying_conditions = all(row1_condition, 1);
    cond_matrix = matrix(:, columns_satisfying_conditions);

    max_damp_bin =  0.01;
    damp_interval = max_damp_bin*0.05;                 % based on 5% change
    damp_bins = 0.001:damp_interval:max_damp_bin;


    [num1,edge2] = histcounts(cond_matrix(2,:), damp_bins);
    [~, indices2] = sort(num1, 'descend');

    damping_Lower_bound = edge2(indices2(1));
    damping_Upper_bound = edge2(indices2(1))+damp_interval;

    row2_condition = cond_matrix(2, :) >= damping_Lower_bound & ...
        cond_matrix(2, :) <= damping_Upper_bound;


    cond_matrix2 = [cond_matrix2,cond_matrix(:,find(all(row2_condition,1),1,'first'))];

end


end


%% A function for Normalize and real(ize) the mode shapes

function Real=ComplexModeToRealMode(Complex)
% This function converts the complex mode shape to the real valued one
% Reference: Operationa modal analysis of civil engineering structures page
% 182 and 183
% Rotate the complex mode shapes (see the RotationMaxCor function, below)
Complex=RotationMaxCor(Complex);
% 1: find the modulus of the mode shape
Modul=abs(Complex);

% 2: normalize the modulus
if nargin < 2
    Modul=Modul/max(Modul);
end

% 3: Find the phase of each component
Phase=angle(Complex);

% 4: find the sign of each component
for I=1:size(Complex,1)
    if Modul(I,1)~=0
        Sign(I,1)=SignFinder(Phase(I));
    else
        Sign(I,1)=0;
    end
end

% 5: compute the real valued mode shape
Real=Modul.*Sign;

end

function Sign=SignFinder(In)
if In>=0 && In<pi/2
    % First quarter
    Sign=+1;
elseif In>=pi/2 && In<=pi
    % Second quarter
    Sign=-1;
elseif In>=-pi && In<-pi/2
    % Third quarter
    Sign=-1;
elseif In>=-pi/2 && In<0
    % Forth quarter
    Sign=+1;
end
end

function out=RotationMaxCor(In)
% This function computes the maximum correlation line and then rotate the
% data with respect to this correlation
X=real(In);
Y=imag(In);
p=polyfit(X,Y,1);% Fit a first order line to the data
Teta=-atan(p(1)); % angle of maximum correlation line
Rot=[cos(Teta)  -sin(Teta) ; sin(Teta)  cos(Teta)]; % Rotation matrix
for I=1:size(In,1);
    N=Rot*[X(I);Y(I)];
    out(I,1)=N(1)+N(2)*1i;
end
end
