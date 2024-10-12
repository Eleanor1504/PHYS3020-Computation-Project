
%PART e


e = 1;     
N = 20;     
k = 1;      
T_start = 1;  
T_end = 3;    
num_temp_steps = 10;  
num_steps = 10000;  
runs = 10;  


T_arr = linspace(T_start, T_end, num_temp_steps);
T_arr = [T_arr, linspace(T_end, T_start, num_temp_steps)];  


sbar_arr = zeros(length(T_arr), runs, num_steps); 
final_spins = zeros(N, N, length(T_arr));  


init_type = 'rand';  


for t_idx = 1:length(T_arr)
    T = T_arr(t_idx);
    
    for run = 1:runs
        
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
                fp = exp(-delU / (k * T));
                if rand < fp
                    ss(i, j) = -ss(i, j);  
                end
            end
            
            
            sbar_arr(t_idx, run, step) = sum(ss(:)) / (N * N);
        end
        
        
        final_spins(:, :, t_idx) = ss;  
    end
end


figure;
subplot(1, 3, 1);
imagesc(final_spins(:, :, 1));  
colormap(gray);  
title('Final Spins at T=1 (Initial)');
colorbar;
axis equal tight;
 
subplot(1, 3, 2);
imagesc(final_spins(:, :, num_temp_steps));  %
colormap(gray);  
title('Final Spins at T=3');
colorbar;
axis equal tight;

subplot(1, 3, 3);
imagesc(final_spins(:, :, end));  )
colormap(gray);  
title('Final Spins at T=1 (Final)');
colorbar;
axis equal tight;

sgtitle('Spin Arrangements at Different Temperatures');

%THIS WORKS !!!!! PART e
