% classdef Musa < handle
%    properties(Access = public,SetObservable, AbortSet)
%       mapping
%       constellation
%       AllConstellationPairs
%    end
%    methods(Access = public)
  function MusaPar = Musa(User,realization,MusaPar)
      K = numel(User);
      Nc = User{1}.Nc ;
      n = User{1}.N;
      M = User{1}.Waveform_Size;
      k = log2(M);
%       receivertype = ReceiverType;
      for user_index_n =1:K
%         = randi([0 1],numberofbits*N,1);
%        User{user_index_n}.ConvertSentDataBitstoInteger();
      
       dataIn = randi([0 1],n*k,1);  % Generate vector of binary dat
       dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
       dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers
       User{user_index_n}.Sent_Data_k = dataSymbolsIn;
       dataModG = qammod(dataSymbolsIn,M); % Gray coding, phase offset = 0  
       ak_original{user_index_n} = dataModG;
       all_const_point = qammod(0:M-1,M);
       alfa = sum(abs(all_const_point).^2)/M;
%        scatterplot(dataModG)
       
       %Spreading Sequence 
       
       if(mod(realization,50)==0 || realization ==1)
       m = 3;
       vect = [-ceil(m/2)+1:ceil(m/2)-1];
       sk_book = vect + j*vect'; 
       sk_book = reshape(sk_book,1,numel(sk_book));
       r = randi([1,m^2],1,Nc);
       for kk=1:Nc
       sk_(kk) = sk_book(r(kk));
       end
%        load('Code4.mat')
%        sk{user_index_n} = bestcode(:,user_index_n);
       sk{user_index_n} = sk_.';
       beta = sum(abs(sk{user_index_n}).^2);
       sk{user_index_n} = (1/(sqrt(beta)*sqrt(alfa))) * sk{user_index_n};
%        sk{user_index_n} = bestcode(:,user_index_n);
       User{user_index_n}.Codebook = all_const_point.*sk{user_index_n};
       scaling_factor =  sqrt(1/(sum(sum(abs(User{user_index_n}.Codebook).^2))/M));
       User{user_index_n}.Codebook = scaling_factor.*User{user_index_n}.Codebook ;
       end
       
      
      end
      
if(mod(realization,50)==0 || realization ==1)
   MusaPar.scodes = sk;
   MusaPar.alfa = alfa;
   MusaPar.constellation = all_const_point;
   MusaPar.ak_original = ak_original; 
end
%        r_ = perms(1:size(V_map,1));
%        Ored = bitor(r_ == 1 , r_ == 2);
%        [~,uniqueindex] = unique(Ored,'rows');
%        r_ = r_(uniqueindex,:);


%        r_ = perms(1:size(V_map,1));
%        sk =          
  
  
  end
      
%    end
%    
% end