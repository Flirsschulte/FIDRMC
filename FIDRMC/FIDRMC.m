%**************************************************************************
% Fuzzy Impulse Noise Detection and Reduction Method for Colour Images
%
%  The paper of the FIDRMC is proposed in: 
%
%  Stefan Schulte, Val�rie De Witte, , Mike Nachtegael
%  Dietrich Van Der Weken and  Etienne E. Kerre:
%  A Fuzzy Two-step Filter for Impulse Noise Reduction From Colour Images,
%  IEEE Transactions on Image Processing 15(11), 2006, 3567 - 3578
%  
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 15/01/06
%
% Inputs:  A = the noisy input image
%          O = the original noise-free image
% Outputs:  Fs = the filtered image 
%**************************************************************************
function Fs = FIDRMC(A,O)
   [M,N,DIM] = size(A);
   A = double(A);
   A1 = A(:,:,1);
   A2 = A(:,:,2);
   A3 = A(:,:,3);
   
   F1 = double(A1);
   F2 = double(A2);
   F3 = double(A3);
   
   flag = 0; loop = 0;
   while (flag == 0) 
       fixed1 = Detection(A1);
       if ((fixed1(1,1) ~= -1) || (loop > 4))
           flag = 1;
       else
           loop = loop + 1;
       end
   end
   flag = 0; loop = 0;
   while (flag == 0) 
       fixed2 = Detection(A2);
       if ((fixed2(1,1) ~= -1) || (loop > 4))
           flag = 1;
       else
           loop = loop + 1;
       end
   end
   flag = 0; loop = 0;
   while (flag == 0) 
       fixed3 = Detection(A3);
       if ((fixed3(1,1) ~= -1) || (loop > 4))
           flag = 1;
       else
           loop = loop + 1;
       end
   end
      
   Windowsize = 1;
   [F1,F2,F3] = DenoiseFIDRMC(A1,A2,A3,fixed1,fixed2,fixed3,Windowsize);
   
   Fs(:,:,1) = F1;   
   Fs(:,:,2) = F2;   
   Fs(:,:,3) = F3;   
   if ((fixed1(1,1)== -1) && (fixed2(1,1)== -1) && (fixed3(1,1)== -1))
      disp(sprintf('======================='));
      Fs = FIDRMCRANDOM(A,O);
      if (MSEC(Fs,O,2) > 800)
          Fs = FIDRMCRANDOM(Fs,O);
      end
      disp(sprintf(' Random-Valued Impulse Noise'));
      disp(sprintf('  => Output MSE:  %10.5f',MSEC(Fs,O,2)));
      disp(sprintf('  => Output PSNR: %10.5f',(log(255^2/MSEC(Fs,O,2))/log(10)*10)));
      disp(sprintf('======================='));
   else
      disp(sprintf('======================='));
      disp(sprintf(' Fixed-Valued Impulse Noise'));
      disp(sprintf('  => Output MSE:  %10.5f',MSEC(Fs,O,2)));
      disp(sprintf('  => Output PSNR: %10.5f',(log(255^2/MSEC(Fs,O,2))/log(10)*10)));
      disp(sprintf('======================='));
   end
   