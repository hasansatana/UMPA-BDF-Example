classdef User_Class_eqsame < matlab.System & matlab.system.mixin.Propagates & matlab.system.mixin.Nondirect
   properties(Access = public,SetObservable, AbortSet)
      % _k at kth time interval
      Waveform_Size
      N = 1 
      G
      V_Matrix 
      Constellation_Mapping
      Sent_Data_k
      Guess_Data
      Codeword_k
      MF_Coeff
      Channel_MF_Coeff
      MF_Output_k
      Channel_Response_k 
      CombinedMFCoeff
      CNoise_Power = 1;
      Receiver_Noise_Power = 1;
      Noise
      Codebook 
      Receiver_Correlated_Noise
      R_mn_ij_
      OneuserReceivedData
      ManualMF 
      Observation
      NOMA_samples
      NOMA_symbols
      TransmitData
      Traditional_Channel_Response
      N_cp = 0;
      Receiver_Correlated_NoiseCMF
      Nantenna
      NumberofActiveMultiPathComponent = 0;
      Nc 
      Receiver_Correlated_NoiseCMF_SIC = [];
      MMSE_SIC_Guess_Data_k;
      %%%
         angular_spread_deltal = 10 % degree;
         steering_vector_u = [];
         ro_coeff = 1; 
 
      %%%
   end
   methods(Access = public)
          function obj = User_Class_eqsame(M,n,Lc,N_antenna)
%           global lc
%           global L
%           global Nc
%           global Tc
          Tc = 1 ;
          Nc = n ;
          lc = Lc;
          L = ceil((lc + Nc -1)/Nc);
          obj.Nantenna = N_antenna;
          obj.Nc = Nc;
          obj.Waveform_Size = M ;
%           addlistener(obj,'Sent_Data_k','PostSet',@obj.handlePropEvents);
          end
          
          function ConvertSentDataBitstoInteger(obj)
            data = reshape(obj.Sent_Data_k,[log2(obj.Waveform_Size),ceil(obj.N)]);
            data = data'; % temporary operation
            data = bi2de(fliplr(data));
            obj.Sent_Data_k = data;
              
          end
         
%           function ArrangeSentData(obj)
%                addlistener(obj,'Sent_Data_k','PostSet',@obj.handlePropEvents);
%           end
          function Make_Receiver_Correlated_Noise(obj,lc,Nc)
%            global lc
%            global Nc

           for s = 1:obj.Nantenna 
           obj.Receiver_Correlated_NoiseCMF(:,s) = conv(obj.Noise(:,s),obj.Channel_MF_Coeff(:,s));
           end
           obj.Receiver_Correlated_NoiseCMF = sum(obj.Receiver_Correlated_NoiseCMF,2);
           obj.Receiver_Correlated_NoiseCMF_SIC = obj.Receiver_Correlated_NoiseCMF(lc:end);
%            end
           for s = 1:size(obj.MF_Coeff,2)
           obj.Receiver_Correlated_Noise(:,s) = conv( obj.Receiver_Correlated_NoiseCMF,obj.MF_Coeff(:,s));
           end
                
           obj.Receiver_Correlated_Noise = obj.Receiver_Correlated_Noise((Nc+lc):end,:);
           % DownSampling is done here
           obj.Receiver_Correlated_Noise = obj.Receiver_Correlated_Noise(1:Nc:end,:);
          end
            function Make_MFilterCoeff(obj)
                signalssquare = abs(obj.Codebook(:,1:end)).^2 ;  % No need for normalisation step, because constellation energy is already 1 
                channelsquare = abs(obj.Channel_Response_k).^2;
                %beta = 1./sqrt(sum(signalssquare));             % Güçleri 1 olmuyor o sebeple s?k?nt? var 
                beta_channel = 1./sqrt(sum(channelsquare));
                %obj.MF_Coeff =  beta.*flipud(conj(obj.Codebook));
                obj.MF_Coeff =  flipud(conj(obj.Codebook));
               % obj.Channel_MF_Coeff = beta_channel.*flipud(conj(obj.Channel_Response_k));
                obj.Channel_MF_Coeff =flipud(conj(obj.Channel_Response_k));
            end
         
            function   Make_Codebook(obj)
                 for m = 1:size(obj.Constellation_Mapping,2)
                 obj.Codebook(:,m) =  obj.V_Matrix*obj.Constellation_Mapping(:,m);
                 end
            end
           
%             function handlePropEvents(obj,src,evnt)
%              switch src.Name
%                 case 'Sent_Data_k'
%                     data = reshape(obj.Sent_Data_k,[log2(obj.Waveform_Size),ceil(obj.N)]);
%                     data = data'; % temporary operation
%                     data = bi2de(fliplr(data));
%                    % data = 2*ones(10,1) % deneme case 
%                     evnt.AffectedObject.Sent_Data_k = data;
%                  case 'MF_Coeff'   % buraya u?ram?yor 
%                      for k = 1:size(obj.Noise,1)
%                          for l = 1:size(obj.Noise,2)
%                      obj.Receiver_Correlated_Noise(k,l) = filter(obj.MF_Coeff(:,k),[1],obj.Noise(k,l));
%                          end
%                      end
%              end
%             end
            
            function ManualMF_outputs(obj)
%                  N_cp = 40;
               TransmitData =  obj.Codebook(:,obj.Sent_Data_k(1:end)+1);
               TransmitData = reshape(TransmitData,numel(TransmitData),1);
               TransmitData_wCP = [TransmitData(end-obj.N_cp+1:end); TransmitData];
%                ReceivedData =  conv(obj.Channel_Response_k,TransmitData)
                ReceivedData =  conv(TransmitData_wCP,obj.Channel_Response_k);  %conv2 de bir ibnelik var çevirmeliyim
%                 ReceivedData = ReceivedData(numel(obj.Channel_Response_k)/2-1:end-numel(obj.Channel_Response_k)/2+1)
%                 ReceivedData =  ReceivedData(numel(obj.Channel_Response_k)/2+1:end-numel(obj.Channel_Response_k)/2+1);
%                 ReceivedData =  conv(obj.Channel_Response_k,TransmitData_wCP)
                obj.OneuserReceivedData = ReceivedData ;  %reshape(ReceivedData,numel(ReceivedData),1);
%                for waveform = 1:obj.Waveform_Size
%                MFilteredData{waveform} =  conv2(ReceivedData(:,1:end),obj.CombinedMFCoeff(:,waveform));
%                [~,maxindex] =  max(abs(MFilteredData{waveform}));
%                MFOutput(waveform,:) = MFilteredData{waveform}(maxindex);  
%                end 
            end
            
            
            function Equalizer(obj,r_MF)
               r_MF = r_MF(obj.Nc+1:end);
               DFT_channel_matrix = diag(fft(obj.Channel_Response_k),numel(r_MF));
               NOMA_fft_symbols = pinv(DFT_channel_matrix).*r_MF ;
               NOMA_samples_wcp = ifft(NOMA_fft_symbols);
               obj.NOMA_samples = NOMA_samples_wcp(obj.N_cp+1:end);
               
               NOMA_samples_seperate = reshape(obj.NOMA_samples,4,numel(obj.NOMA_samples)/4);

               for k = 1:length(NOMA_samples_seperate)
                  for i = 1:16
                   Distance(k,i) = norm(NOMA_samples_seperate(:,k)- obj.Codebook(:,i));
                  end
                  [~,symindex] = min(Distance(k,:));
                  Sembol(k) = symindex -1 ;
              end
               
            end
            
   end
         events 
              DataAssigned
         end
   end
