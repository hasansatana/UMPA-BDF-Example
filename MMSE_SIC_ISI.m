function biterrorrate = MMSE_SIC_ISI(User,MusaPar)
% 108 ini de bulmayaa çal???yor ?u an ona göre 
Numberofuser = size(User,2) ;
Lc = size(User{1}.Channel_Response_k,1);
T = numel(User{1}.Sent_Data_k);
Nc = numel(MusaPar.scodes{1});
N0 = 1/1000;
L = ceil((Lc + Nc -1)/Nc);
M = User{1}.Waveform_Size;

for i = 1:Numberofuser
    SentDataTotal{i} = [];
    for k=1:T
    SentData_Symbol{i,k} = User{i}.Codebook(:,User{i}.Sent_Data_k(k)+1);
    SentDataTotal{i} = [SentDataTotal{i} ; SentData_Symbol{i,k}];
    end
end


    r_v = zeros(Numberofuser,Numberofuser);
    for i=1:Numberofuser
        for j=1:Numberofuser
            for N = 1: size(User{i}.Channel_Response_k,2)
                if(Lc==1)
                aa = xcorr(User{i}.Channel_Response_k(:,N),User{j}.Channel_Response_k(:,N),0);
                else
                aa = xcorr(User{i}.Channel_Response_k(:,N),User{j}.Channel_Response_k(:,N)); 
                end
                r_autocorr{i,j}(:,N) = aa;
            end
            r_autocorr{i,j} =  sum(r_autocorr{i,j},2);
             r_h{i,j} = r_autocorr{i,j}(Lc);

             r_v(j,i) = r_h{i,j};
             
        end
      
    end
user_vector = 1:Numberofuser;

 y  = zeros(Lc+numel(SentDataTotal{1})-1,N);
  for user_index_v = 1:Numberofuser
        for N_antenna = 1:N
        y_vs{user_index_v}(:,N_antenna) = conv(SentDataTotal{user_index_v},User{user_index_v}.Channel_Response_k(:,N_antenna));
        end
        y = y  +  y_vs{user_index_v} ; % +vk s? var 
  end

for user_index_u = 1:Numberofuser
    for N_antenna = 1:N
        y_hN{user_index_u}(:,N_antenna) = conv(User{user_index_u}.Channel_MF_Coeff(:,N_antenna),y(:,N_antenna));
    end
        y_h{user_index_u} = sum(y_hN{user_index_u},2) ;
        y_h{user_index_u} = y_h{user_index_u}(Lc:end-Lc+1) + User{user_index_u}.Receiver_Correlated_NoiseCMF_SIC;
end

    
    for user_index_u = 1: Numberofuser     
            for k = 1:T
             y_u{k,user_index_u} = y_h{user_index_u}(Nc*(k-1)+1:Nc*k);   %noise eklendi, kanaldan geçti
            end
    end
    
    counter = 0 ;
    for user_index_v =1:Numberofuser
        g{user_index_v} = [];
        for index = 1:Numberofuser
                g{user_index_v} = [g{user_index_v};r_v(index,user_index_v).*MusaPar.scodes{user_index_v}];
        end
         counter = counter +1 ;
         G(:,user_index_v) = g{user_index_v};
    end
   

for iteration=1:Numberofuser 

     if(iteration>1)
    founduser_index = find(user_vector ==founduser);
    user_vector = [user_vector(1:founduser_index-1) user_vector(founduser_index+1:end)];
     end
    for k = 1:T
        r{k}= [];
        for user_index_u = 1: Numberofuser        
             if iteration>1
                 y_subtract{k,founduser} = r_h{founduser,user_index_u}*MusaPar.scodes{founduser}* a(k,founduser);
                 y_u{k,user_index_u} =  y_u{k,user_index_u} - y_subtract{k,founduser};
             end
         r{k} =  [r{k};y_u{k,user_index_u}];
        end
    end
    Es = MusaPar.alfa;

    for index = 1:Numberofuser

          Correlations_Matrix{index} = zeros(Numberofuser-1);
          SINR(index) = -1000;
       
    end
 
    GW = zeros(size(G));
    GW(:,user_vector) = G(:,user_vector);
     Correlations_MatrixW= GW'*GW; 
     Correlations_MatrixW = Correlations_MatrixW*Es + N0*eye(Numberofuser);
%     
     W_G = GW*Correlations_MatrixW^-1*Es;
         
%%%%%%%%%%%%%%%%% SINR Calculation

    for user_index_m =user_vector

        GSNR = zeros(size(G)); 
       
        GSNR(:,user_vector) = G(:,user_vector);
        GSNR(:,user_index_m) = zeros(Nc*Numberofuser,1);
       
        GSNR(:,find(all(GSNR==0))) = [];
       Correlations_Matrix{user_index_m} = GSNR*GSNR'; 
       Correlations_Matrix{user_index_m} = Correlations_Matrix{user_index_m}*Es + N0*eye(Nc*Numberofuser);
       
       SINR(user_index_m) =  (Es*abs(W_G(:,user_index_m)'*g{user_index_m}).^2)/(abs(W_G(:,user_index_m)'*Correlations_Matrix{user_index_m}*W_G(:,user_index_m)));
       SINR(user_index_m) =  10*log10(SINR(user_index_m));

    end
    [maxSNR, maxSNR_index] = max(SINR);

     founduser = maxSNR_index;

    for k=1:T 
%         a(k,founduser) = W{founduser}' * r{k}/Es;  % Fazladan bi Es bölmesi istiyor nerden bulamad?m 
        a_(k,:) = W_G'*r{k};
        a(k,founduser) =  a_(k,founduser);
    end
   
end


for user_index_m = 1:Numberofuser
      for k=1:T 
       [~,minindex] = min(abs(MusaPar.constellation -  a(k,user_index_m)));
       User{user_index_m}.MMSE_SIC_Guess_Data_k(k) = minindex -1  ;
      
      end
end

numberofbiterror = zeros(1,Numberofuser);
for i = 1:Numberofuser
    
error = find((User{i}.Sent_Data_k(L:end-L+1) - User{i}.MMSE_SIC_Guess_Data_k(L:end-L+1)) ~=0);   %%%%%% ------------------ BURASI DE???MEL? 
fprintf('User Hata %d = %d\n',i,numel(error));
biterror = fliplr(de2bi(User{i}.Sent_Data_k(L:end-L+1),log2(M))) - fliplr(de2bi(User{i}.MMSE_SIC_Guess_Data_k(L:end-L+1),log2(M)));
biterror = reshape(biterror,numel(biterror),1);
numberofbiterror(i) = numel(find(biterror ~=0));
%  end
totalbiterror = sum(numberofbiterror);
N = T;
biterrorrate = totalbiterror/(log2(M)*(N-2*(L-1))*Numberofuser);




end