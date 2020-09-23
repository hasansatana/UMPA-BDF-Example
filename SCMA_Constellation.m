classdef SCMA_Constellation < handle
   properties(Access = public,SetObservable, AbortSet)
      mapping
      constellation
      AllConstellationPairs
   end
   methods(Access = public)
      function obj = SCMA_Constellation(constellation)
         obj.constellation = constellation;
         obj.SelectConstellation();
         global hQAMMod_RE1;
         global hQAMMod_RE2;
         hQAMMod_RE1 = comm.GeneralQAMModulator; 
         hQAMMod_RE2 = comm.GeneralQAMModulator; 
         hQAMMod_RE1.Constellation = obj.mapping(1,:);
         hQAMMod_RE2.Constellation = obj.mapping(2,:);
         scatterplot(hQAMMod_RE1.Constellation);
         xlabel('Re'); ylabel('Im');  title('Selected Constellation Mapping RE1');  grid on;  text(real(hQAMMod_RE1.Constellation)+0.01, imag(hQAMMod_RE1.Constellation), dec2bin(0:numel(hQAMMod_RE1.Constellation)-1,log2(numel(hQAMMod_RE1.Constellation))));
         scatterplot(hQAMMod_RE2.Constellation);
         xlabel('Re'); ylabel('Im'); title('Selected Constellation Mapping RE2') ; grid on; text(real(hQAMMod_RE1.Constellation)+0.01, imag(hQAMMod_RE1.Constellation), dec2bin(0:numel(hQAMMod_RE2.Constellation)-1,log2(numel(hQAMMod_RE2.Constellation)))); 
         obj.AllConstellationPairs = [ hQAMMod_RE1.Constellation(1:end) ; hQAMMod_RE2.Constellation(1:end) ];
      end
      
     function   SelectConstellation(obj)
        switch obj.constellation
        case 1  %LDS 
         obj.mapping = [0.5+0.5j 0.5-0.5j -0.5+0.5j -0.5-0.5j; 0.5+0.5j 0.5-0.5j -0.5+0.5j -0.5-0.5j];
        case 2 %T4QAM
         obj.mapping = [3/sqrt(10)+0j  -1/sqrt(10)+0j 1/sqrt(10)+0j -3/sqrt(10)+0j;1/sqrt(10)+0j 3/sqrt(10)+0j -3/sqrt(10)+0j -1/sqrt(10)+0j];
        case 3  %4LQAM
         obj.mapping = [-sqrt(2)/2  -sqrt(2)/2 sqrt(2)/2 sqrt(2)/2;-sqrt(2)*i/2 -sqrt(2)*i/2 sqrt(2)*i/2 sqrt(2)*i/2];
        case 4 %4BAO
         obj.mapping = [0.5019-0.4981j  +0.5019+0.4981j -0.5019-0.4981j -0.5019+0.4981j; 0.5019-0.4981j -0.5019-0.4981j 0.5019+0.4981j -0.5019+0.4981j];
        case 5 %4CQM
          obj.mapping = [1  0 0 -1; 0 1 -1 0];
        case 6 %4BEKO
          obj.mapping = [-0.7586-0.1274i  , 0.2626- 0.8822i , -0.2583 + 0.6244i, 0.7543+0.3852j;-0.1835+0.6120i, 0.1707 - 0.3517i ,0.4121-0.6113i, -0.3993+ 0.3509j];
        case 7  %16 LDS 
          obj.mapping = [-0.75-0.75j -0.75-0.25j -0.75+0.75j -0.75+0.25j -0.25-0.75j -0.25+0.75j -0.25-0.25j -0.25+0.25j  0.75-0.75j 0.75-0.25j 0.75+0.75j 0.75+0.25j 0.25-0.75j 0.25-0.25j 0.25+0.75j 0.25+0.25j];
          obj.mapping = [obj.mapping; obj.mapping];
        case 8
           obj.mapping = [ 3*sqrt(5)/10+3j*sqrt(5)/10  3*sqrt(5)/10-sqrt(5)*j/10  3*sqrt(5)/10+sqrt(5)*j/10  3*sqrt(5)/10-3j*sqrt(5)/10 -sqrt(5)/10+3j*sqrt(5)/10 -sqrt(5)/10-j*sqrt(5)/10 -sqrt(5)/10+j*sqrt(5)/10 -sqrt(5)/10-3j*sqrt(5)/10 ...
             sqrt(5)/10+3j*sqrt(5)/10 sqrt(5)/10-j*sqrt(5)/10 sqrt(5)/10+j*sqrt(5)/10 sqrt(5)/10-3j*sqrt(5)/10 -3*sqrt(5)/10+3j*sqrt(5)/10 -3*sqrt(5)/10-j*sqrt(5)/10 -3*sqrt(5)/10+j*sqrt(5)/10 -3*sqrt(5)/10-3j*sqrt(5)/10 ...
             ;  sqrt(5)/10+j*sqrt(5)/10  sqrt(5)/10+3j*sqrt(5)/10 sqrt(5)/10-3j*sqrt(5)/10 sqrt(5)/10-j*sqrt(5)/10  3*sqrt(5)/10+j*sqrt(5)/10 3*sqrt(5)/10+3j*sqrt(5)/10 3*sqrt(5)/10-3j*sqrt(5)/10 3*sqrt(5)/10-j*sqrt(5)/10 -sqrt(5)/10+j*sqrt(5)/10 -3*sqrt(5)/10+j*sqrt(5)/10 -3*sqrt(5)/10+3j*sqrt(5)/10 -3*sqrt(5)/10-3j*sqrt(5)/10 -3*sqrt(5)/10-j*sqrt(5)/10   -sqrt(5)/10+3j*sqrt(5)/10  -sqrt(5)/10-3j*sqrt(5)/10 -sqrt(5)/10-j*sqrt(5)/10];
        case 9   %16LQAM 
           obj.mapping = [ sqrt(2)/2+j*sqrt(2)/2  sqrt(2)/2 sqrt(2)/2-j*sqrt(2)/2 j*sqrt(2)/2  0 0 sqrt(2)/2 -j*sqrt(2)/2 j*sqrt(2)/2 0 0  -j*sqrt(2)/2 -sqrt(2)/2+j*sqrt(2)/2  -sqrt(2)/2 -sqrt(2)/2  -sqrt(2)/2-j*sqrt(2)/2 ...
             ;0 j*sqrt(2)/2  -j*sqrt(2)/2 0 sqrt(2)/2 sqrt(2)/2+j*sqrt(2)/2  sqrt(2)/2-j*sqrt(2)/2  sqrt(2)/2 -sqrt(2)/2 -sqrt(2)/2+j*sqrt(2)/2 -sqrt(2)/2-j*sqrt(2)/2 -sqrt(2)/2  0  j*sqrt(2)/2 -j*sqrt(2)/2 0 ]
         end
    end

   end
end