function [parameters] = readAcqDataInfo(fname)
%This function opens a text file saved by Eric's version of the labview
%program and gets the scan parameters that were set by the user before
%clicking acquire. 
%The list of what each number in paramter means can be displayed by
%uncommenting line 34. This is the output of that line.
%     {'vend'             }
%     {'vstart'           }
%     {'hend'             }
%     {'hstart'           }
%     {'#vPixel'          }
%     {'#hPixel'          }
%     {'Origin Position X'}
%     {'Origin Position Y'}
%     {'End Position X'   }
%     {'End Position Y'   }
%     {'EMExposure'       }
%     {'ZPos'             }
%     {'XCount'           }
%     {'YCount'           }
%     {'ZCount'           }
%     {'XDis'             }
%     {'YDis'             }
%     {'ZDis'             }
%     {'FSR'              }
%     {'GHz/pixel'        }
%     {'Laser Y'          }
%     {'Laser X'          }
%   author: mnikolic@umd.edu
%%
% fname='C:\Users\mniko\Documents\MATLAB\data\20181120 3D aligned matrix\a2\11_20_04_04_45PM AcqData.tdms';
fname=[fname(1:end-5),'Info.txt'];
T=readtable(fname,'ReadRowNames',true);
%T.Row
% parameters=T.Variables;
parameters=T;
end

