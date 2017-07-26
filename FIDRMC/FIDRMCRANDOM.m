%**************************************************************************
% Fuzzy Impulse Noise Detection and Reduction Method for Colour Images
%
%  The paper of the FIDRMC is proposed in: 
%
%  Stefan Schulte, Valérie De Witte, , Mike Nachtegael
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



function Fs = FIDRMCRANDOM(A,O)
  [M,N,DIM] = size(A);
   A = double(A);
   A1 = A(:,:,1);
   A2 = A(:,:,2);
   A3 = A(:,:,3);
   
   FR = double(A1);
   FG = double(A2);
   FB = double(A3);

   a = 100.0; b = 190.0;
   Windowsize = 1;
   F = double(A);
   F2 = double(A);

   begin = 0;
   loop = 1;
   while (begin == 0)
      for k = 0:4 
         F = double(F2);
         FR = double(F(:,:,1));
         FG = double(F(:,:,2));
         FB = double(F(:,:,3));
         mem(M,N) = 0;
         
         memR= Detectionb(FR,a/(2^k),b/(2^k));
         memG= Detectionb(FG,a/(2^k),b/(2^k));
         memB= Detectionb(FB,a/(2^k),b/(2^k));
         [FR2,FG2,FB2] = DenoiseFIDRMC2(FR,FG,FB,memR,memG,memB,Windowsize);
          F2(:,:,1) = FR2;
          F2(:,:,2) = FG2;
          F2(:,:,3) = FB2;

         if (MSEC(F,O,3) < MSEC(F2,O,3)) | (loop == 40)
            begin = 1;
            break;
         end
      end
      loop = loop +1;
   end
   Fs = F;   
