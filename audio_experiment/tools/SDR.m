function out = SDR(yd,y,allow_delay)
%   out = SDR(y_d,y,allow_delay);
%   Inputs:
%           yd  - Distorted signal vector length x 1 | 1 x length
%           y   - Original not-distorted signal vector length x 1 | 1 x length
%           allow_delay - true | false to allow delay into computations 
%   Outputs:
%           out - Signal-to-Distortion ratio in dB

    if(ndims(yd)~=2 || numel(yd) ~= length(yd) || ndims(y)~=2 || numel(y) ~= length(y) || ...
             length(y) ~= length(yd))
        error('Unexpected input dimensions!');
    end
    
    if(size(yd,2)==1)
        yd = yd.';
    end
    
    if(size(y,2)==1)
        y = y.';
    end
    
    if(allow_delay)
        [yd,y,D] = alignsignals(yd,y,[],'truncate');
        if(abs(D)>=(length(y)/10))
            warning('Allowed delay seems to be larger than 1/10 length of input signals.');
        end
    end  
    
    alpha = (yd*y')/(y*y');

    out = (sum(yd.*conj(yd))/sum((yd-(y.*alpha)).*conj(yd-(y.*alpha))));
end
