function output=calcAreas(F,PSD,VLF,LF,HF,flagNorm)
%calcAreas - Calulates areas/energy under the PSD curve within the freq
%bands defined by VLF, LF, and HF. Returns areas/energies as ms^2,
%percentage, and normalized units. Also returns LF/HF ratio.
%
%Inputs:
%   PSD: PSD vector
%   F: Freq vector
%   VLF, LF, HF: array containing VLF, LF, and HF freq limits
%   flagNormalize: option to normalize PSD to max(PSD)
%Output:
%
%Usage:
%   
%
%   Modified from Gary Clifford's ECG Toolbox: calc_lfhf.m   
    if nargin<6
       flagNorm=false;
    end
    
    %normalize PSD if needed
    if flagNorm
        PSD=PSD/max(PSD);
    end
    % find the indexes corresponding to the VLF, LF, and HF bands
    iVLF= (F>=VLF(1)) & (F<=VLF(2));
    iLF = (F>LF(1)) & (F<=LF(2));
    iHF = (F>HF(1)) & (F<=HF(2));
      
    % calculate raw areas (power under curve), within the freq bands (ms^2)
    %aVLF=trapz(F(iVLF),PSD(iVLF));
    %aLF=trapz(F(iLF),PSD(iLF));
    %aHF=trapz(F(iHF),PSD(iHF));
    aVLF=trapz(F.*iVLF,PSD.*iVLF);
    aLF=trapz(F.*iLF,PSD.*iLF);
    aHF=trapz(F.*iHF,PSD.*iHF);
    %Michele Castiotta
    %Provo a mettere il * 1000000 in quanto mi trovo sfasato di scala
    aVLF = aVLF * 1000000; %单位ms*ms/Hz
    aLF = aLF * 1000000;
    aHF = aHF * 1000000;
    
    aTotal=aVLF+aLF+aHF;
    
    %calculate normalized areas (relative to HF+LF, n.u.)
    nLF=aLF/(aLF+aHF);
    nHF=aHF/(aLF+aHF);
    
    %Michele Castriotta
    %moltiplico * 100 
    nLF = nLF * 100;
    nHF = nHF * 100;
    
    %calculate LF/HF ratio
    lfhf =aLF/aHF;
            
    %create output structure
    if flagNorm
        output.aVLF=(aVLF*1000)/1000;
        output.aLF=(aLF*1000)/1000;
        output.aHF=(aHF*1000)/1000;
        output.aTotal=(aTotal*1000)/1000;
    else
        output.aVLF=(aVLF*100)/100; % round
        output.aLF=(aLF*100)/100;
        output.aHF=(aHF*100)/100;
        output.aTotal=(aTotal*100)/100;
    end
    
    output.pVLF=(pVLF*10)/10;
    output.pLF=(pLF*10)/10;
    output.pHF=(pHF*10)/10;
    output.nLF=(nLF*1000)/1000;
    output.nHF=(nHF*1000)/1000;
    output.LFHF=(lfhf*1000)/1000;
end