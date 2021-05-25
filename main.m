clear all      
close all
% Author: Hasan Aykut Satana / ASELSAN / Middle East Technical University
% Algorithm: Dr. Gökhan Muzaffer Güvensen - Hasan Aykut Satana 
%##############################Introduction##########################################
% This is the main BER/AIR calculation code for an examplary MIMO NOMA
% uplink transmission in long dispersive channels with single carrier modulations. 
% You can select the receiver type MMSE-SIC or UMPA-BDF or combination of them.
% You can also arrange a receiver to obtain single user matched filter performance of the receiver for a given uplink NOMA scenerio
% Code domain NOMA schemes which are SCMA and MUSA can be selected for the
% simulations. 
% AIR calculations are based on Generalized Mututal Information concept
% For the simulations along with the BER/AIR with SNR,  BER/AIR vs Loading curves can also be obtained
%where loading is number of users over code length 
% (GMI) 
% Receivers are given in .mex format. whose copyrights are belong to
% authors

%#############################Parameters############################################
% N : Number of antennas 
% K : Number of users 
% Nc: Code length
% Lc: Channel length 
% M : Modulation order
%###################################################################################

% Receiver Types
%1) UMPA-BDF 2) MMSE-SIC 3) MMSE-SIC aided UMPA-BDF 4) UMPA-BDF with Matched Filter Bounds  5) UMPA-BDF with Multi-User Bounds
% "UMPA"  "MMSESIC"  "UMPAwMMSESIC" "UMPAMFB"  "UMPAMU"
ReceiverType = 'UMPA';
OperationMode = 'BERvSNR'; % BERvSNR or BERvLoading
MusaParameters = struct();
if(strcmp(OperationMode,'BERvLoading'))
    MinNumberofUsers = 1;
    MaxNumberofUsers = 32; 
    UserIncrement = 2;
    UserVector = MinNumberofUsers:UserIncrement:MaxNumberofUsers;
elseif(strcmp(OperationMode,'BERvSNR'))
    NumberofUsers = 16; 
    UserVector = NumberofUsers;
else
    printf('No valid selection can be found.');
end 
%################################ Variable Simulation Parameters ##########
    NOMACodeSelection = 'MUSA';            % SCMA or MUSA 
    NumberofRealization = 1;               % #of realization in each SNR 
    NumberofRealizationRefreshCodes = 50;  % Change usercodes 
    NAntenna =4;                           % #of antenna in the receiver 
    M = 4;                                 % Alphabet (Constellation) size
    Nc = 8;                                % Codeword size, Codelength
    dm = 2;                                % SCMA Dimension, # of nonzero elements    
    Lc = 16 ;                              % Channel length
    N = 100;                               % Number of NOMA symbols is one frame 2* Lc / Nc  + 1 
%#####################################################################################################################

for LoadingIndex = 1:numel(UserVector) % if loading=false then only 1 loop 
    
    K = UserVector(LoadingIndex) ;         % Number of users according to loading scenerio K/Nc
    k = 0 ;                                %Start time index
    numberofbits = log2(M);                % # of bits in each symbol
    L = ceil((Lc + Nc -1)/Nc);             % Effective channel length    
    N = N + 2*(L-1);                       % #of NOMA symbols in one packet together with tranining symbols added 
    SNRCounter = 0;                        % SNR counter
    User = cell(K,1);                      % cell for User objects
    SNRList = 0:2:30;                     % SNR List
        
    for user_index_n=1:K
        User{user_index_n} = User_Class(M,Nc,Lc,NAntenna);
        User{user_index_n}.N = N;
    end


    for SNRdB = SNRList
        SNRCounter = 1+ SNRCounter;
        SNR = 10^(SNRdB/10);                   % Linear value of SNR
        N0 = 1/SNR;                            % Noise power
        %N0 = 0;
        BitErrorRate = [];
        biterrorrateMMSE = [];
        
        for Realization = 1:NumberofRealization
            fprintf('SNRdB = %d dB,  Realization = %d\n',SNRdB,Realization);
            
            n_k = wgn(N*Nc,NAntenna,N0,'linear','complex');    % Uncorrelated upsampled receiver noise
            r_MF = 0;                                          % MF outputs initialization
            
            if(strcmp(NOMACodeSelection,'MUSA'))
                MusaParameters = Musa(User,Realization,MusaParameters);
            elseif(strcmp(NOMACodeSelection,'SCMA'))
                %SCMA()
            end

            for user_index_n=1:K
            %   rng('shuffle')
                power = sqrt(1/NAntenna/Lc);                   %Noise power per antenna*channel tab
                User{user_index_n}.Channel_Response_k = power*wgn(Lc,NAntenna,1,'linear','complex');
                User{user_index_n}.Make_MFilterCoeff();
                User{user_index_n}.Noise = n_k;
                User{user_index_n}.Receiver_Correlated_Noise = [];
                User{user_index_n}.Make_Receiver_Correlated_Noise(Lc,Nc); % Correlated Noise after matched filters for each user
                r_MF = User{user_index_n}.OneuserReceivedData + r_MF;
            end
            
            kstart = + L  ;                   % Since first L is for training
            kend = N  - L + 1 ;               % Last L is also for training as we have 1 packet 

            for user_index_n=1:K
                User{user_index_n}.Guess_Data = zeros(N,1); 
                User{user_index_n}.Guess_Data(1:kstart-1) = User{user_index_n}.Sent_Data_k(1:kstart-1); % Training symbol placement
                User{user_index_n}.Guess_Data(kend+1:end) = User{user_index_n}.Sent_Data_k(kend+1:end); % Training symbol placement
                User{user_index_n}.MMSE_SIC_Guess_Data_k = zeros(N,1);
            end
            % "UMPA"  "MMSESIC"  "UMPAwMMSESIC" "UMPAMFB"  "UMPAMU"
            if(strcmp(ReceiverType,'MMSESIC')||strcmp(ReceiverType,'UMPAwMMSESIC'))
%                 ARGS_MMSE{1} = coder.typeof(User,[],true);
%                 ARGS_MMSE{2} = coder.typeof(MusaParameters,[],true);
                biterrorrateMMSE_ = MMSE_SIC(User,MusaParameters);
                biterrorrateMMSE(Realization) = biterrorrateMMSE_;
            end
            
 %################################ Define User Structure
 %######################## Since 'codegen' is not working with class arrays

            for user_index_n=1:K
                UserStr = struct();
                UserStr.Codebook =  User{user_index_n}.Codebook;
                UserStr.Waveform_Size = M; 
                UserStr.Lc = Lc;
                UserStr.Nc = Nc;
                UserStr.N = N - 2*(L-1);
                UserStr.L = L;
                UserStr.k_start = kstart;
                UserStr.k_end = kend;
                UserStr.N0 = N0;
                UserStr.Guess_Data = User{user_index_n}.Guess_Data;
                UserStr.Sent_Data_k = User{user_index_n}.Sent_Data_k;
                UserStr.MMSE_SIC_Guess_Data_k = User{user_index_n}.MMSE_SIC_Guess_Data_k;
                UserStrCellArray{user_index_n} = UserStr;     
            end
           %########################Trick for code generation to make
           %variable size base ##########################################s
                UserStrCellArray_type = UserStrCellArray;
                UserStrCellArray_type{K+1} =  UserStr;
                UserStrCellArray_type{K+1}.Codebook = complex(zeros(Nc*6,M*8));
                UserStrCellArray_type{K+1}.Guess_Data = zeros(N+200,1);
                UserStrCellArray_type{K+1}.Sent_Data_k = zeros(N+200,1);
                UserStrCellArray_type{K+1}.MMSE_SIC_Guess_Data_k = zeros(N+200,1);
               
%##########################################################################################################
           
% "UMPA"  "MMSESIC"  "UMPAwMMSESIC" "UMPAMFB"  "UMPAMU"
            if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
                if(Realization == 1 || mod(Realization,NumberofRealizationRefreshCodes)==0)
                  
                    ARGS = coder.typeof(UserStrCellArray_type,[inf,inf],true);
%                     codegen CalculateUserCodeCorrelations -args ARGS
                    tic
                    fprintf('Correlations are being calculated...\n');
                    Gmnij = CalculateUserCodeCorrelations_mex(UserStrCellArray);
                    fprintf('Correlations are ready...\n');
                    toc
                end
%############################################################################################################
                %####################### This is where the user
                %correlations with MIMO NOMA codes are calculated. ######
                for user_m = 1:K
                    for user_n = 1:K
                         for l_count = 1:2*L-1  % l_count is the channel tap 
                                for waveform_i = 1:M
                                    for waveform_j = 1:M
                                        Rmnij{user_m,user_n,l_count}(waveform_i,waveform_j) = trace(User{user_n}.Channel_Response_k'*Gmnij{user_m,user_n,waveform_i,waveform_j,l_count}*User{user_m}.Channel_Response_k);          
                                    end
                                end       
                        end
                    end
                end

                 for user = 1:K
                        for k=1:N
                                c  = zeros(M,1);  % c indicates which NOMA code is selected in the User codebook 
                                c(User{user}.Sent_Data_k(k)+1) = 1;
                                ckm{user,k} = c; % for all users and for all the symbols in the packet c variable is created
                        end
                 end
%##################################################################################################################################
%################################# This is where the MF outputs are
%calculated for all users and for all symbol time.
                for k= kstart:kend  %k indicates the symbol time .
                    for user_n = 1:K
                        r = 0 ;
                         for user_m = 1:K
                              l_count = 1 ;
                              for l= -L +1 : L-1 
                                      r = r +  Rmnij{user_m,user_n,l_count}.'*ckm{user_m,k- l};
                                      l_count = l_count + 1;
                              end
                         end
                         rnMF{user_n,k} = r ;  % MF output for user n
                         rnMFnoisy{user_n,k} = r + User{user_n}.Receiver_Correlated_Noise(kstart,:).'; % Correlated Noise added (due to MFs) 
                         rnMF{user_n,k} =    rnMFnoisy{user_n,k};
                    end
                end
 %############################################ MF outputs are ready now ###############################################
 %############################################## Decoding Part###########        
                Rmnij_fake = Rmnij;
                Rmnij_fake{K,K,2*L} = zeros(64*M,64*M);
 
                 if(Realization == 1 && SNRCounter ==1)
                    ARGS_2{1} = coder.typeof(UserStrCellArray_type,[inf,inf],true);
                    ARGS_2{2} = coder.typeof(Rmnij_fake,[],true);
                    ARGS_2{3} = coder.typeof(rnMF,[],true);
%                     codegen ReceiverIterations -args ARGS_2
                 end

                 if(strcmp(ReceiverType,'UMPA'))
                    [biterrorrateit,GMIit] = ReceiverIterations_mex(UserStrCellArray,Rmnij,rnMF);
                elseif(strcmp(ReceiverType,'UMPAwMMSESIC'))
                    [biterrorrateit,GMIit] = ReceiverIterations_with_MMSESIC_mex(UserStrCellArray,Rmnij,rnMF);
%                 elseif(strcmp(ReceiverType,'UMPAMFB'))
%                     [biterrorrate_it,GMIit] =  ReceiverIterationsSU_mex(UserStrCellArray,Rmnij,rnMF);
%                 elseif(strcmp(ReceiverType,'UMPAMU'))
%                     [biterrorrate_it,GMIit] =  ReceiverIterationsMU_mex(UserStrCellArray,Rmnij,rnMF);
                 else
                     printf('There is no such option, neither it is not implemented nor invalid')
                 end
                BitErrorRate(:,Realization) = biterrorrateit;
                GMI_realization(:,Realization) = GMIit;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        biterrorrateSNR(:,SNRCounter) = mean(BitErrorRate,2);
        biterrorrateMMSESNR(SNRCounter) = mean(biterrorrateMMSE);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        biterrorrateOverloading(LoadingIndex) = biterrorrateSNR(end,SNRCounter);
        if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')||strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
        GMISNR(:,SNRCounter) = mean(GMI_realization,2)/Nc;
        GMIOverloading(LoadingIndex) = GMISNR(end,SNRCounter) ;
        end
    end
end
%############################### Plots-- BER/GMI vs SNR/Overloading ##########################################

%###################First Plots in terms of SNR##########
iteration=6;
if(strcmp(OperationMode,'BERvSNR'))
    if(strcmp(ReceiverType,'UMPAMU')) iteration= 5; end
    if(strcmp(ReceiverType,'UMPAMFB')) iteration= 1; end
%######3 GMI- SNR Plots  (Maximum uncoded data rate) #################
    figure;
    for i =1:iteration
        plot(SNRList,GMISNR(i,:))
        title(['GMI-SNR-M',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(NAntenna)]) 
        xlabel('E_s/N_0 (dB)')
        ylabel('Average Uncoded Data Rate per User')
        grid minor
        hold on
    end
    saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'fig.fig'])
    saveas(gcf,['GMI_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'tif.tif'])

savestring = ['Data',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'.mat'];
save(savestring);
%#####################3 BER / SNR Plots #################################
figure
hold on
if(strcmp(ReceiverType,'MMSESIC')||strcmp(ReceiverType,'UMPAwMMSESIC'))
semilogy(SNRList,biterrorrateMMSESNR)
end

if(strcmp(ReceiverType,'UMPA')||strcmp(ReceiverType,'UMPAwMMSESIC')|| strcmp(ReceiverType,'UMPAMFB')||strcmp(ReceiverType,'UMPAMU')) 
semilogy(SNRList,biterrorrateSNR)
end

xlabel('E_s/N_0 (dB)')
ylabel('BER')
title(['BER',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(NAntenna)])
grid minor 

end
%#####################################################################################33
%###################Second Plots in terms of Overloading##########
%################  L O A D I N G ############ 
if(strcmp(OperationMode,'BERvLoading'))
    Overloading_vector = UserVector./Nc;  
%#################################GMI vs Overloading %%%%%%%%%%%%%%%%%%%%%
    figure;
    plot(Overloading_vector,GMIOverloading)
    title(['GMI-SNR-M',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(NAntenna),'.mat']) 
    grid minor
    saveas(gcf,['GMI_Overloading_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'fig.fig'])
    saveas(gcf,['GMI_Overloading_SNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'tif.tif'])
    savestring = ['OverloadingData',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'.mat'];
    save(savestring);
    
    figure;
    semilogy(Overloading_vector,biterrorrateOverloading)
    xlabel('Es/N0 (dB)')
    ylabel('BER')
    title(['BER',num2str(M),'Lc',num2str(Lc),'N-antenna',num2str(NAntenna)])
    grid minor 
    saveas(gcf,['FigSNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'.fig'])
    saveas(gcf,['TiffSNR_M',num2str(M),'Lc',num2str(Lc),'N_antenna',num2str(NAntenna),'.tif'])
    %############################ #####################################
 end





   
