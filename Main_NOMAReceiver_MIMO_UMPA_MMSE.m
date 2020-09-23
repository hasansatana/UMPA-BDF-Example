clear all      
close all
% tic 
% UMPA  MMSESIC  UMPAwMMSESIC  UMPAwMMSESICmax UMPAMFB  UMPAMU
ReceiverType = 'UMPA'; % 'UMPA'
Loading = 'off';
%Spatial Channel Parameters 
angular_spread_deltal = 10; % degree;
steering_vector_u = [];
ro_coeff = 1; 
NumberofActiveMultiPathComponent = 6;
%%%
MusaPar = struct();
K_vector = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]*2;
K_vector = [20]; 
for Overloading = 1:numel(K_vector)
numberofrealization = 1;
N_antenna =4;
M = 4;          % Alphabet size
Nc = 8;          % Codeword size
dm = 2;          % Dimension 
K = K_vector(Overloading) ;     % Number of users
Lc = 32 ; 
T = Nc; 
Tc = 1;          % 1 symbol time 
k = 0 ;  %start time
numberofbits = log2(M);

L = ceil((Lc + Nc -1)/Nc);
N = 100;  % 2* Lc / Nc  + 1 
N = N + 2*(L-1);
SNR_count = 0;
for user_index_n=1:K
User{user_index_n} = User_Class_eqsame(M,Nc,Lc,N_antenna);
User{user_index_n}.N = N;
end


SNRGeneral = 30:2:30;
for SNR_dB = SNRGeneral

SNR_count = 1+ SNR_count;
SNRdB = SNR_dB; % dB
SNR = 10^(SNRdB/10);
N0 = 1/SNR;
% N0 = 0;
biterrorrate = [];
biterrorrate_MMSE_real = [];
for realization = 1:numberofrealization


fprintf('SNRdB = %d,  Realization = %d\n',SNRdB,realization);

n_k = wgn(N*Nc,N_antenna,N0,'linear','complex'); 

%Nc + Lc - 1 +
r_MF = 0;

MusaPar = Musa(User,realization,MusaPar);

for user_index_n=1:K
%   rng('shuffle')
 power = sqrt(1/N_antenna/Lc);
  User{user_index_n}.Channel_Response_k = power*wgn(Lc,N_antenna,1,'linear','complex');
  User{user_index_n}.Make_MFilterCoeff();
  User{user_index_n}.Noise = n_k;
  User{user_index_n}.Receiver_Correlated_Noise = [];
  User{user_index_n}.Make_Receiver_Correlated_Noise(Lc,Nc);  %& açýlacak
  r_MF = User{user_index_n}.OneuserReceivedData + r_MF;
  %%%%%%%%%
end

k_start = + L  ; 
k_end = N  - L + 1 ;

for user_index_n=1:K
User{user_index_n}.Guess_Data = zeros(N,1); 
User{user_index_n}.Guess_Data(1:k_start-1) = User{user_index_n}.Sent_Data_k(1:k_start-1);
User{user_index_n}.Guess_Data(k_end+1:end) = User{user_index_n}.Sent_Data_k(k_end+1:end);
User{user_index_n}.MMSE_SIC_Guess_Data_k = zeros(N,1);
end
% UMPA  MMSESIC  UMPAwMMSESIC  UMPAwMMSESICmax UMPAMFB  UMPAMU
if(strcmp(ReceiverType,'MMSESIC')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax'))
biterrorrate_MMSE = MMSE_SIC_ISI(User,MusaPar);
biterrorrate_MMSE_real(realization) = biterrorrate_MMSE;
end

for user_index_n=1:K
Userstruct = struct;
Userstruct.Codebook = User{user_index_n}.Codebook;
Userstruct.Waveform_Size = M; 
Userstruct.Lc = Lc;
Userstruct.Nc = Nc;
Userstruct.N = N - 2*(L-1);
Userstruct.L = L;
Userstruct.k_start = k_start;
Userstruct.k_end = k_end;
Userstruct.N0 = N0;
Userstruct.Guess_Data = User{user_index_n}.Guess_Data;
Userstruct.Sent_Data_k = User{user_index_n}.Sent_Data_k;
Userstruct.MMSE_SIC_Guess_Data_k = User{user_index_n}.MMSE_SIC_Guess_Data_k;
UserCellArray{user_index_n} = Userstruct;
end

% UMPA  MMSESIC  UMPAwMMSESIC  UMPAwMMSESICmax UMPAMFB  UMPAMU
if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
lc = Lc;
Waveform_Size = M;
ck = [zeros(Waveform_Size,L-1,K) zeros(Waveform_Size,N,K) zeros(Waveform_Size,L-1,K)];
user_m =1 ;
R_ij = {};
     
if(realization == 1 || mod(realization,50)==0)
fprintf('Correlations are calculated\n');
% ARGS = coder.typeof(UserCellArray);
% codegen CalculateCorrelation -args ARGS
tic
G_mn_ij = CalculateCorrelation_mex(UserCellArray);
toc
end

for user_m = 1:K
    for user_n = 1:K
                 for l_count = 1:2*L-1
                        for waveform_i = 1:M
                            for waveform_j = 1:M
                             Rmnij{user_m,user_n,l_count}(waveform_i,waveform_j) = trace(User{user_n}.Channel_Response_k'*G_mn_ij{user_m,user_n,waveform_i,waveform_j,l_count}*User{user_m}.Channel_Response_k);          
                            end
                        end
%                   R_mn_ij{,user_n,:,waveform_j,l_count} = Rmnij(:,:);            
                end
    end
end

 for user = 1:K
        for k=1:N
                c  = zeros(M,1);
                c(User{user}.Sent_Data_k(k)+1) = 1;
                ck_m{user,k} = c;
        end
 end

for k= k_start:k_end
    for user_n = 1:K
            r = 0 ;
             for user_m = 1:K
                      l_count = 1 ;
                      for l= -L +1 : L-1 
                              r = r +  Rmnij{user_m,user_n,l_count}.'*ck_m{user_m,k- l};
                              l_count = l_count + 1;
                      end
             end
             r_nMF{user_n,k} = r ;
             r_nMF_noisy{user_n,k} = r + User{user_n}.Receiver_Correlated_Noise(k_start,:).';
             r_nMF{user_n,k} =    r_nMF_noisy{user_n,k};
    end
end
%    

T_indicator = 1 ;
P_Ik = 1/M;
vk_nm = cell(K,K,k_end);
uk_nm = cell(K,K,k_end);
for i=1:numel(vk_nm)
    vk_nm{i} = (zeros(M,1));  %n için   %m th user n th user  k th time 
end
for i=1:numel(uk_nm)
    uk_nm{i} = (zeros(M,1));  %n için   %m th user n th user  k th time 
end

% vksum = zeros(M,1);
% l0 = ceil(Lc/4+1 - 1/4);
     %%%%%%%%%%%%%%% End of  Time Domain Forward-Backward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 


%  if(realization == 1 && SNR_count ==1)
%    t = coder.typeof(r_nMF);
%    t = makeHomogeneous(t);
%    ARGS_2{1} = coder.typeof(UserCellArray);
%    ARGS_2{2} = coder.typeof(Rmnij);
%    ARGS_2{3} =  t; %coder.typeof(r_nMF); 
% %    cfg =coder.config('lib');
% % UMPA  MMSESIC  UMPAwMMSESIC  UMPAwMMSESICmax UMPAMFB  UMPAMU
%     if(strcmp(ReceiverType,'UMPA'))
%         codegen ReceiverIterations -args ARGS_2 
%         codegen ReceiverIterationsnofor -args ARGS_2 
%     elseif(strcmp(ReceiverType,'UMPAwMMSESIC'))
%         codegen  ReceiverIterations_with_MMSESIC -args ARGS_2 
%     elseif(strcmp(ReceiverType,'UMPAwMMSESICmax'))
%         codegen ReceiverIterations_with_MMSESIC_max -args ARGS_2 
%     elseif(strcmp(ReceiverType,'UMPAMFB'))
%          codegen ReceiverIterationsSU -args ARGS_2
%     elseif(strcmp(ReceiverType,'UMPAMU'))
%         codegen ReceiverIterationsMU -args ARGS_2
%     end
   
%  end
%  tic
 if(strcmp(ReceiverType,'UMPA'))
        [biterrorrate_it,GMI_it] = ReceiverIterations_mex(UserCellArray,Rmnij,r_nMF);
%         [biterrorrate_it2,GMI_it2] = ReceiverIterationsnofor_mex(UserCellArray,Rmnij,r_nMF);
    elseif(strcmp(ReceiverType,'UMPAwMMSESIC'))
        [biterrorrate_it,GMI_it] = ReceiverIterations_with_MMSESIC_mex(UserCellArray,Rmnij,r_nMF);
    elseif(strcmp(ReceiverType,'UMPAwMMSESICmax'))
        [biterrorrate_it,GMI_it] = ReceiverIterations_with_MMSESIC_max_mex(UserCellArray,Rmnij,r_nMF);
    elseif(strcmp(ReceiverType,'UMPAMFB'))
        [biterrorrate_it,GMI_it] =  ReceiverIterationsSU_mex(UserCellArray,Rmnij,r_nMF);
    elseif(strcmp(ReceiverType,'UMPAMU'))
        [biterrorrate_it,GMI_it] =  ReceiverIterationsMU_mex(UserCellArray,Rmnij,r_nMF);
 end
 
%   [biterrorrate_it,GMI_it] =  ReceiverIterations(UserCellArray,Rmnij,r_nMF);
%     [biterrorrate_it,GMI_it] = ReceiverIterations_SoftC(UserCellArray,Rmnij,r_nMF);
%  toc

   %%
% end
% if(strcmp(ReceiverType,'UMPA'))
biterrorrate(:,realization) = biterrorrate_it;
GMI_realization(:,realization) = GMI_it;

% end
end
end

%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

biterrorrateSNR(:,SNR_count) = mean(biterrorrate,2);
% biterrorrateSNR2(:,SNR_count) = mean(biterrorrate2,2);
biterrorrate_MMSE_SNR(SNR_count) = mean(biterrorrate_MMSE_real);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
biterrorrateOverloading(Overloading) = biterrorrateSNR(end,SNR_count);
if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
GMI_SNR(:,SNR_count) = mean(GMI_realization,2)/Nc;
GMI_Overloading(Overloading) = GMI_SNR(end,SNR_count) ;
end
end


end

iteration = 6;
if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
if(strcmp(ReceiverType,'UMPAMFB'))
    iteration = 1 ;
elseif(strcmp(ReceiverType,'UMPAMU'))
    iteration = 5;
else 
    iteration = 6;
end

% if(strcmp(Loading,'off'))
figure;
for i =1:iteration
plot(SNRGeneral,GMI_SNR(i,:))
title(['GMI-SNR-M',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(N_antenna),'.mat']) 
grid minor
hold on
end
saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'fig.fig'])
saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'tif.tif'])
end
% Overloading_vector = K_vector./8;
% %%%
string = ['Data',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'.mat'];
save(string);
figure
if(strcmp(ReceiverType,'MMSESIC')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax'))
semilogy(SNRGeneral,biterrorrate_MMSE_SNR)
end
hold on
if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
semilogy(SNRGeneral,biterrorrateSNR)
end

% figure;
% if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAwMMSESICmax')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
% semilogy(SNRGeneral,biterrorrateSNR2)
% end

xlabel('Es/N0')
% ylabel('BER')
title(['BER',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(N_antenna),'.mat'])
grid minor 
% for i = 1:iteration
% legendtitle{i} = (['iteration-',num2str(i)]);
% end
% legend(legendtitle);

if(strcmp(Loading,'on'))
%################  L O A D I N G ############ 

iteration = 6;
figure;
Overloading_vector = K_vector./Nc;



plot(Overloading_vector,GMI_Overloading)
title(['GMI-SNR-M',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(N_antenna),'.mat']) 
grid minor
saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'fig.fig'])
saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'tif.tif'])

string = ['Data',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'.mat'];
save(string);

Overloading_vector = K_vector./Nc;
%%%
string = ['Data',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'.mat'];
save(string);
figure
semilogy(Overloading_vector,biterrorrateOverloading)
xlabel('Es/N0')
ylabel('BER')
title(['BER',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(N_antenna),'.mat'])
grid minor 
for i = 1:iteration
legendtitle{i} = (['iteration-',num2str(i)]);
end
%############################ #####################################
end
saveas(gcf,['FigSNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'.fig'])
saveas(gcf,['TiffSNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(N_antenna),'.tif'])




   
