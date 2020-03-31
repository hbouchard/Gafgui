function    T = fct_SignaltoTransmission(S,LUT,channel)

err = 0;
if channel==1
    T = interp1(double(LUT.raw),double(LUT.red),double(S))/(2^16-1);
elseif channel==2
    T = interp1(double(LUT.raw),double(LUT.green),double(S))/(2^16-1);
elseif channel==3
    T = interp1(double(LUT.raw),double(LUT.blue),double(S))/(2^16-1);
else
   err = 1; 
end