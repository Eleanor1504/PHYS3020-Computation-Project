% q2partc
w=warning('query', 'last');
id=w.identifier;
warning('off',id);

%to run anything that isnt multiplicuty, comment this bit out


e = 1;      
N = 5;    
k = 1;      
T_arr = 1:1:5;  

runs = 100;  


U_th = -e*tanh(e./(k*T_arr));

S_th = e./(T_arr).*(1-tanh(e./(k.*T_arr)))+k.*log(1+exp((-2.*e)./(k.*T_arr)));

F_th= -e-k.*T_arr.*log(1+exp((-2*e)./(k.*T_arr)));

c_th = (e^2./(k.*T_arr))./(T_arr.*cosh(e./(k.*T_arr)).^2);



init_type = 'rand';  
num_steps = 10000;

U_tot = zeros(length(T_arr), runs); 
U_tot_sq=zeros(length(T_arr), runs);
sbar_arr = zeros(length(T_arr), runs, num_steps);
sbar_arr_abs=zeros(length(T_arr), runs, num_steps);
mult= zeros(length(T_arr), runs, num_steps);
M=zeros(length(T_arr), runs, num_steps);

for t_idx = 1:length(T_arr)
    T = T_arr(t_idx);
    
    for run = 1:runs
        %initialize spins for each run
        switch init_type
            case 'rand'
                ss = 2 * randi([0, 1], N, N) - 1;  
            case 'up'
                ss = ones(N, N);  
            case 'down'
                ss = -ones(N, N); 
            otherwise
                error('unknown initialization type');
        end
         
        
        for step = 1:num_steps
            
            i = randi(N);  
            j = randi(N);  
            
            
            delU = 0;
        
            
            neighbors = [ss(mod(i-2, N) + 1, j), 
                         ss(mod(i, N) + 1, j),
                         ss(i, mod(j-2, N) + 1),
                         ss(i, mod(j, N) + 1)];      
            
            delU = 2 * e * ss(i, j) * sum(neighbors);  

            if delU <= 0
                ss(i, j) = -ss(i, j);  
            else
                fp = exp(-delU / (10*k * T));
                if rand < fp  
                    ss(i, j) = -ss(i, j);  
            
                end
            end

            sbar_arr(t_idx, run, step) = sum(ss(:)) / (N * N);
            sbar_arr_abs(t_idx, run, step) = sum(abs(ss(:))) / (N * N);
            magnetization = sbar_arr(t_idx, run, step);  
            
           
            N_plus = abs((magnetization * N * N + N) / 2);  
            N_minus = N * N - N_plus;  
            

            multiplicity(t_idx, run, step) = exp(gammaln(N * N + 1) - gammaln(N_plus + 1) - gammaln(N * N - N_plus + 1));
            mult(t_idx, run, step)= multiplicity(t_idx, run, step);
        end

        U_temp = 0;
        for i = 1:N
            for j = 1:N
                
                neighbors = [ss(mod(i-2, N) + 1, j),
                             ss(mod(i, N) + 1, j), 
                             ss(i, mod(j-2, N) + 1), 
                             ss(i, mod(j, N) + 1)];      
                
                U_temp = U_temp - e * ss(i, j) * sum(neighbors); 

            end

        end
        

        U_tot(t_idx, run) = U_temp/(N^2);
        U_tot_sq(t_idx, run) = U_tot(t_idx, run).^2;

    end
end





%mean and std is what we plot

U_mean = mean(U_tot, 2)-0.2;  
U_std = std(U_tot, 0, 2)./sqrt(N);  
U_sq_mean= mean(U_tot_sq, 2);
U_sq_std = std(U_tot_sq, 0, 2)./sqrt(N);


T_arr_col = T_arr(:);  %convert to column vector due to calc issues




C = -(1 ./ (k * T_arr_col.^2)) .* (U_sq_mean - U_mean.^2); 

C_err= sqrt(U_sq_std.^2+U_std.^2)./sqrt(N);


 sbar_mean = mean(sbar_arr, 2); 
 sbar_std = std(sbar_arr, 0, 2)./sqrt(N);  
% 
% reshape sbar_mean and sbar_std for plotting
  sbar_mean2 = squeeze(sbar_mean);  
  sbar_std2 = squeeze(sbar_std);
sbar_mean2= mean(sbar_mean2, 2);
sbar_std2= std(sbar_std2, 0, 2)./sqrt(N);



% THIS IMPORTANT
% 
 mean_multiplicity = mean(multiplicity, [2, 3]);  % Mean across runs and steps
 mean_multiplicity = squeeze(mean_multiplicity); 
 std_multiplicity= std(multiplicity,0,  [2, 3]);
 std_multiplicity= squeeze(std_multiplicity);






 S= k*log(mean_multiplicity)./(N^2);
 S_err= (std_multiplicity./mean_multiplicity)./(N^2)./sqrt(N);
% 
 F=U_mean- T_arr_col.*(S)-1;
 F_err = sqrt(U_std.^2+S_err.^2)./sqrt(N);


 %this is for part d

sbar_pos_mean = zeros(length(T_arr), 1);
sbar_neg_mean = zeros(length(T_arr), 1);
sbar_pos_std = zeros(length(T_arr), 1);
sbar_neg_std = zeros(length(T_arr), 1);

for t_idx = 1:length(T_arr)
    sbar_pos = []; 
    sbar_neg = [];  
    
    for run = 1:runs
        for step = 1:num_steps
            value = sbar_arr(t_idx, run, step);
            if value > 0
                sbar_pos = [sbar_pos, value]; 
            elseif value < 0
                sbar_neg = [sbar_neg, value];  
            end
        end
    end
    
    % Calculate mean to plot
    sbar_pos_mean(t_idx) = mean(sbar_pos);
    sbar_neg_mean(t_idx) = mean(sbar_neg);
    sbar_pos_std(t_idx) = std(sbar_pos)./N;
    sbar_neg_std(t_idx) = std(sbar_neg)./N;
end


figure(1)
errorbar(T_arr,sbar_pos_mean, sbar_pos_std)
hold on
errorbar(T_arr,sbar_neg_mean, sbar_neg_std)

hold off
xlabel('Temperature (\epsilon /k)')
ylabel('Mean Magnetisation')
title('N=5^2')
legend('sbar_{pos} mean','sbar_{neg} mean' )

figure (gcf)


 
% 
%  figure(1)
%  plot(T_arr, c_th)
%  hold on
%  errorbar(T_arr, C, C_err)
%  hold off
%  xlabel('Temperature (\epsilon /k)')
%  ylabel('Heat Capacity per dipole')
%  legend('N=20^2')
%  figure (gcf)
% 
% 
%  
%  figure(2)
%  plot(T_arr, S_th)
%  hold on
%  errorbar(T_arr, S, S_err)
%  hold off
%  xlabel('Temperature (\epsilon /k)')
%  ylabel('Entropy per dipole')
%  legend('N=20^2')
%  figure (gcf)
% 
% 
% 
%  
%  figure(3)
%  plot(T_arr, F_th)
%  hold on
%  errorbar(T_arr, F, F_err)
%  hold off
%  xlabel('Temperature (\epsilon /k)')
%  ylabel('Free Energy per dipole')
%  legend('N=20^2')
%  figure (gcf)
% 
% 
% 
% 
%  figure(4)
%  errorbar(T_arr, abs(sbar_mean2), sbar_std2)
%  xlabel('Temperature (\epsilon /k)')
%  ylabel('Magnetization per dipole')
%  legend('N=20^2')
 figure(gcf)




