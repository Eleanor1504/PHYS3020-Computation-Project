
e=1;

N= 500;

% k=1.38*10^(-23);
k=1;

runs = 200;
% T=2;

% T_arr= [0.01, 1, 2, 6];
% T_arr = [0.001:0.01:20];

% T_arr = [0.5, 1, 2];
T=2;

% T_arr= [1:10:1000];

U_th = -e*tanh(e./(k*T_arr));


init_type = 'rand'; %'rand', 'up', down';


switch init_type
    case 'rand'
        ss = 2*randi([0,1], 1, N) -1;
     
    case 'up'
        ss = ones(size(N));

    case 'down'
        ss = -1.*ones(size(N));
    otherwise
        error('unknown')

end


% U = -e*
% 
% E_0

s = randsample(ss, 1);
% 
% sumspin = circshift(spinchain, [0,1])+circshift(spinchain, [0, -1]);
% 
% Em= -e.*spinchain.*sumspin;
% Em;
% E=0.5.*sum(Em)
% 
    % Emean= E/length(spinchain)
ss0=ss;
Uarray=zeros(1,length(T_arr));

delU=zeros(1,N);
sbararr= zeros(1, length(T_arr));
% shistarr=zeros(1, 1000)
% shisarr = zeros (1, runs);
% 
%    for l=1:length(runs)
%         for j=1:10000
%     %     shist= shistarr(l);    
%             i = randi(N);
%         %     for i=1:N
%                 if (i==1)
%                     delU(i) = -e*0.5* (ss(N)*ss(i) +ss(i)*ss(i+1));
%                 elseif (i==N)
%                      delU(i) = -e*0.5* (ss(i-1)*ss(i) +ss(i)*ss(1));
%                 else
%                      delU(i) = -e*0.5* (ss(i-1)*ss(i) +ss(i)*ss(i+1));
%                  end
%         %     end
%             if delU(i)<=0
%                 ss(i)=-ss(i);
%             else
%                 fp = exp(-10*delU(i)/(k*T));
%                 if rand>fp
%                     ss(i)=-ss(i);
%                 end
%             end
%         end
%     shisarr(l)= sum(ss)/N;
%     end 

shisarr = zeros(1, runs);

% Monte Carlo simulation
for l = 1:runs
    for j = 1:100000  % Number of Monte Carlo steps
        i = randi(N);  % Randomly select a spin to flip
        
        % Calculate change in energy (delU)
        if i == 1
            delU = -e * 0.5 * (ss(N) * ss(i) + ss(i) * ss(i + 1));
        elseif i == N
            delU = -e * 0.5 * (ss(i - 1) * ss(i) + ss(i) * ss(1));
        else
            delU = -e * 0.5 * (ss(i - 1) * ss(i) + ss(i) * ss(i + 1));
        end

        % Metropolis acceptance criterion
        if delU <= 0
            ss(i) = -ss(i);  % Accept the flip
        else
            fp = exp(-delU / (k * T_arr(1)));  % Use T_arr(1) for this example
            if rand > fp
                ss(i) = -ss(i);  % Accept the flip with probability
            end
        end
    end
    
    % Store the average magnetization for this run
    shisarr(l) = sum(ss) / N;
end 
