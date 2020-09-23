function biterrorrate = MMSE_SIC(User,MusaPar)
% 108 ini de bulmayaa çal???yor ?u an ona göre 
Numberofuser = size(User,2) ;
Lc = size(User{1}.Channel_Response_k,1);
T = numel(User{1}.Sent_Data_k);
Nc = numel(MusaPar.scodes{1});
N0 = 1/1000;
L = ceil((Lc + Nc -1)/Nc);
M = User{1}.Waveform_Size;
for i = 1:Numberofuser
    for k=1:T
    SentData_Symbol{i,k} = User{i}.Codebook(:,User{i}.Sent_Data_k(k)+1);
    end
end


    r_v = zeros(Numberofuser,Numberofuser);
    for i=1:Numberofuser
        for j=1:Numberofuser
            for N = 1: size(User{i}.Channel_Response_k,2)
                aa = xcorr(User{i}.Channel_Response_k(:,N),User{j}.Channel_Response_k(:,N),0);
                r_autocorr{i,j}(:,N) = aa;
            end
            r_autocorr{i,j} =  sum(r_autocorr{i,j},2);
             r_h{i,j} = r_autocorr{i,j}(Lc);

             r_v(i,j) = r_h{i,j};
        end

    end
user_vector = 1:Numberofuser;

    for k = 1:T
        r{k}= [];
        for user_index_u = 1: Numberofuser
            v_u{user_index_u} = User{user_index_u}.Receiver_Correlated_NoiseCMF_SIC;
            v_uk{k,user_index_u} = v_u{user_index_u}((k-1)*Nc+1:k*Nc);
            y_u{k,user_index_u} = zeros(Nc,1);
            for user_index_v = 1:Numberofuser
                y_vs{k,user_index_u,user_index_v} = r_h{user_index_u,user_index_v}*SentData_Symbol{user_index_v,k};
                y_u{k,user_index_u} = y_u{k,user_index_u}  +  y_vs{k,user_index_u,user_index_v} ; % +vk s? var 
            end
             y_u{k,user_index_u} =  y_u{k,user_index_u} +  v_uk{k,user_index_u};   
        end
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
                 y_subtract{k,founduser} = r_h{user_index_u,founduser}*MusaPar.scodes{founduser}* a(k,founduser);
                 y_u{k,user_index_u} =  y_u{k,user_index_u} - y_subtract{k,founduser};
             end
         r{k} =  [r{k};y_u{k,user_index_u}];
        end
    end

    % g_v  = zeros(Nc*Numberofuser,Numberofuser);

    counter = 0 ;

    for user_index_v =1:Numberofuser
        g{user_index_v} = [];
        for index = 1:Numberofuser
                g{user_index_v} = [g{user_index_v};r_v(index,user_index_v).*MusaPar.scodes{user_index_v}];
        end
         counter = counter +1 ;
    end

    Es = MusaPar.alfa;

    for index = 1:Numberofuser
          Correlations{index} = zeros(Nc*Numberofuser); 
          SNR(index) = -1000;
    end

    for user_index_m =user_vector

        for i = user_vector
            if i ~= user_index_m
            Correlations{user_index_m} = Correlations{user_index_m} +  g{i}*g{i}';
            end
        end
       Correlations{user_index_m} =Correlations{user_index_m}*Es + N0*eye(Nc*Numberofuser)*N0;
       SNR(user_index_m) = g{user_index_m}'*Correlations{user_index_m}^-1*g{user_index_m}*Es^2;
       SNR(user_index_m) = 10*log10(real(SNR(user_index_m)));
       
    end
    [maxSNR, maxSNR_index] = max(SNR);

    % W = zeros(Numberofuser*Nc,1);
     founduser = maxSNR_index;
    Correlations_w{founduser} = zeros(Nc*Numberofuser);
    for i = user_vector
       Correlations_w{founduser} = Correlations_w{founduser} + g{i}*g{i}';
    end
     Correlations_w{founduser} =Correlations_w{founduser}*Es + eye(Nc*Numberofuser)*N0;
     W{founduser} = Correlations_w{founduser}^-1*g{founduser}*Es^2;


    for k=1:T 
        a(k,founduser) = W{founduser}' * r{k}/Es;  % Fazladan bi Es bölmesi istiyor nerden bulamad?m 
    end
   
end


for user_index_m = 1:Numberofuser
      for k=1:T 
       [~,minindex] = min(abs(MusaPar.constellation -  a(k,user_index_m)));
       User{user_index_m}.Guess_Data(k) = minindex -1  ;
      
      end
end
% 
% L laz?m M laz?m 
numberofbiterror = zeros(1,Numberofuser);
for i = 1:Numberofuser
error = find((User{i}.Sent_Data_k(L:end-L+1) - User{i}.Guess_Data(L:end-L+1)) ~=0);   %%%%%% ------------------ BURASI DE???MEL? 
fprintf('User %d Hata %d\n',i,numel(error)) 
biterror = fliplr(de2bi(User{i}.Sent_Data_k(L:end-L+1),log2(M))) - fliplr(de2bi(User{i}.Guess_Data(L:end-L+1),log2(M)));
biterror = reshape(biterror,numel(biterror),1);
numberofbiterror(i) = numel(find(biterror ~=0));
%  end
totalbiterror = sum(numberofbiterror);
N = T;
biterrorrate = totalbiterror/(log2(M)*(N-2*(L-1))*Numberofuser);




end