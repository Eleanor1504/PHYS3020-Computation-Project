%part c plots 

e=1;

N= 100;


runs = 10;
% k=1.38*10^(-23);
k=1;

% T=2;

% T_arr= [0.01, 1, 2, 6];
T_arr = [0.5:1:10];

% T_arr = [0.5, 1, 2];

% T_arr= [1:10:1000];

U_th = -e*tanh(e./(k*T_arr));

S_th = e./(T_arr).*(1-tanh(e./(k.*T_arr)))+k.*log(1+exp((-2.*e)./(k.*T_arr)));

F_th= -e-k.*T_arr.*log(1+exp((-2*e)./(k.*T_arr)));

c_th = (e^2./(k.*T_arr))./(T_arr.*cosh(e./(k.*T_arr)).^2);


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

M= zeros (size(sbar_mean));
mult= zeros(size(sbar_mean));

delU=zeros(1,N);
sbararr= zeros(1, length(T_arr));
U_arr=zeros(1,length(T_arr));

U2_arr=zeros(1,length(T_arr));
% Uerr = zeros(10000, length(T_arr));
% sbarerr = zeros(10000, length(T_arr));
% shistarr=zeros(1, 1000)
 
    for k=1:length(T_arr)
        T = T_arr(k);
        for l = 1:runs   
        sbarvals = zeros(100000, 1);
        Uvals= zeros(100000, 1);
        U2vals= zeros (100000, 1);
            for j=1:100000
        %     shist= shistarr(l);    
                i = randi(N);
            %     for i=1:N
                    if (i==1)
                        delU(i) = -e*0.5* (ss(N)*ss(i) +ss(i)*ss(i+1));
                    elseif (i==N)
                         delU(i) = -e*0.5* (ss(i-1)*ss(i) +ss(i)*ss(1));
                    else
                         delU(i) = -e*0.5* (ss(i-1)*ss(i) +ss(i)*ss(i+1));
                     end
            %     end
                if delU(i)<=0
                    ss(i)=-ss(i);
                else
                    fp = exp(-delU(i)/(k*T));
                    if rand>fp
                        ss(i)=-ss(i);
                    end
                end
            sbarvals(j) = sum(ss)/N;
            Uvals(j)=sum(delU-1)/N;
            U2vals(j)= Uvals(j).^2;

%             U_sq_vals(j) = Uvals(j).^2;
            end
        
        sbararr (l, k)= mean(abs(sbarvals));   
        U_arr (l, k)= mean(Uvals);
%         U2_arr(l, k) = mean(Uvals.^2);
        U2_arr(l,k)=mean(U2vals);
        M(l, k)=10*N*sbararr(l,k);
        K(l,k)=round((M(l,k)+N)/2,0);
        mult(l,k)= nchoosek(N, K(l,k));
        end
%         Uarray(k)= sum(delU-1)/N;
%         sbararr(k)= sum(ss)/N;
    end
sbar_mean = mean(sbararr, 1);
U_mean= mean(U_arr, 1);
% U_sq_mean = mean(U_sq_arr, 1);
% U_mean_sq = U_mean.^2;
U2_mean = mean(U2_arr, 1);
mult_mean= mean(mult, 1);


sbar_err=std(sbararr, 0, 1);
U_err= std(U_arr, 0, 1);
U2_err= std(U2_arr, 0, 1);
C_err= sqrt(U2_err.^2+U_err.^2);
mult_err= std(mult, 0, 1);


S = k*log(mult_mean)/(10*N);

S_err= mult_err./mult_mean;

F_err = sqrt(S_err.^2+U_err.^2);

% 
%  for k=1:length(T_arr)
%  
%      M(k)= 10*N*sbar_mean(k);
%      K(k)=round((M(k)+N)/2,0);
%      mult(k)= nchoosek(N, K(k));
%  end


% S = k*log(mult)/(10*N);

% S_err = 

% F= zeros (size(S));
F= U_mean-(T_arr.*S);

C = 20.*(U2_mean-U_mean.^2)./((T_arr).^2);

% for k = 1:length(T_arr)
%     M = round(N * sbar_mean(k));            % Calculate M from average magnetization
%     n_plus = (N + M) / 2;                   % Number of +1 spins
%     % Ensure n_plus is an integer
%     if mod(N + M, 2) ~= 0 || n_plus < 0 || n_plus > N
%         multiplicity(k) = 0;  % Invalid state
%     else
%         multiplicity(k) = nchoosek(N, n_plus); % Calculate multiplicity
%     end
% end
% sbar_err = ones(1,length(sbararr))*std(sbararr);
% U_err = ones(1,length(Uarray))*std(Uarray)


%basic eqns for plots, figure out how to loop over temps and store 10 vals
% Utot=sum(delU);
% 
% u=Utot/N;
