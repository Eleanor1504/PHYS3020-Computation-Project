%q2 part a and b
e = 1;      
N = 100;    
k = 1;      
T =4;      

%initialize spins
init_type = 'rand';  

switch init_type
    case 'rand'
        ss = 2 * randi([0, 1], N, N) - 1;  %change initialisation for 2D
    case 'up'
        ss = ones(N, N); 
    case 'down'
        ss = -ones(N, N); 
    otherwise
        error('unknown');
end

ss0 = ss;
num_steps = 10000;
for step = 1:num_steps
    
    i = randi(N);  
    j = randi(N); 
    
    
    delU = 0;

    % change the periodic boundary conditions for 2d
    %this decides which dipoles neighbour the selected one
    neighbors = [ss(mod(i-2, N) + 1, j),
                 ss(mod(i, N) + 1, j), 
                 ss(i, mod(j-2, N) + 1), 
                 ss(i, mod(j, N) + 1)];      
    

    delU = 2 * e * ss(i, j) * sum(neighbors);


    
    if delU <= 0
        ss(i, j) = -ss(i, j);  
    else %probability fucntion is the same
        fp = exp(-delU / (k * T));
        if rand < fp  
            ss(i, j) = -ss(i, j);  
        end
    end
end


U_tot = 0;

% Sum over spins and neighbour iteractions to find energy contributions
for i = 1:N
    for j = 1:N
        neighbors = [ss(mod(i-2, N) + 1, j),
                     ss(mod(i, N) + 1, j), 
                     ss(i, mod(j-2, N) + 1), 
                     ss(i, mod(j, N) + 1)];      % Right
        
        U_tot = U_tot - e * ss(i, j) * sum(neighbors);  %update total energy
    end
end

% Average energy per spin
u = U_tot / (N * N);


%plot results for each temperature

tiledlayout(1,2)


nexttile
heatmap(ss0)
title('Initla (T=2.0 \epsilon /k)')
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

nexttile
heatmap(ss)
title('Final(T=2.0 \epsilon /k)')
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

figure(gcf)

