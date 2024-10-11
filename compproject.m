%comp project
% syms e
e=1;

N= 100;

k=1;

T=0.5;

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

 s = randsample(ss, 1);


%save initial array
ss0=ss;

%create array of the right size to store each energy contribution
delU=zeros(1,N);


%j tells matlab how many times to run the simultion
%ie how many chances for each dipole to flip


for j=1:10000
    i = randi(N);
%criteria for the periodic boundary conditoins to simulate infinite lattice
        if (i==1)
            delU(i) = -e* (ss(N)*ss(i) +ss(i)*ss(i+1));
        elseif (i==N)
             delU(i) = -e* (ss(i-1)*ss(i) +ss(i)*ss(1));
        else
             delU(i) = -e* (ss(i-1)*ss(i) +ss(i)*ss(i+1));
         end
%flip if dipoles are different
    if delU(i)<=0
        ss(i)=-ss(i);
    else
        fp = exp(-delU(i)/(k*T)); %flip probability, T dependent 
        if rand>fp
            ss(i)=-ss(i); %flip based on the probability
        end
    end
end



%basic eqns for plots, figure out how to loop over temps and store 10 vals


Utot=sum(delU);


u=Utot/N;


%plotting initial and final states
tiledlayout(1,2)


nexttile
heatmap(ss0)
title('Initial State (0.5 \epsilon / k)')
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

nexttile
heatmap(ss)
title('Final State (0.5 \epsilon / k)')
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));



figure(gcf)


