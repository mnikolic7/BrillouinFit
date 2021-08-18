function [ shifts ] = sep2shift( shifts, FSR)
%This function takes in your shift results (from fitting) and converts them
%to the actual Brillouin shift. WARNING: your results must be already
%converted to GHz if you want this to work. 
%   Separation=distance between antistokes to stokes of the previous
%   diffraction order. 
%   shift = Actual Brillouin shift - distance from Raileigh center peak.
%   author: mnikolic@umd.edu

% disp('Input to this function must be in GHz, if you want it to work');
warning('Input to this function must be in GHz, if you want it to work');
shifts=(FSR-shifts)/2;

end

