classdef Gridness<handle
    %GRIDNESS With this class one can analyze the gridcells
    %   TODO
    
    properties
        spikingData         % Contains the spiking data of one cell in the format [timestamp,x,y]
        spikingCoordinates  % Contains the spiking data of one cell in the format [x,y]
    end
    
    methods
        
        function readData(obj,path)
            % This function receives a path to a binary file that stores
            % data in the format timestamp,x,y (float32) . It saves this
            % data to the variable spikingData. It is assumed that the
            % maximum/minimum values in either x or y direction is 1
            % respectivly -1.
            
            fid = fopen(path,'r');
            tmp = fread(fid,'*float32');
            fclose(fid);
            obj.spikingData = [tmp(1:3:end) tmp(2:3:end) tmp(3:3:end)];
            obj.spikingCoordinates = obj.spikingData(:,2:3);
        end
        
        function [rateMap] = calcRateMap(obj,scale)
            % This functions calculates a rate map in the given scale using
            % the spikingCoordinates data. It is assumed that the 
            % maximum/minimum values in either x or y direction is 1 
            % respectivly -1.
            
            % use only positive coordinates, scale and round them so that
            % they are  indexable
            data = (obj.spikingCoordinates +1 )*scale;
            data = round(data);

            % initialize rateMap
            rateMap = zeros(2*scale,2*scale);

            % fill rateMap according to the times the cell fires in that
            % field
            for i=1:size(data)
                 j = data(i,1);
                 k = data(i,2);
                 rateMap(j,k) = rateMap(j,k) + 1;
            end
        end
        
        function [autoCorrelationMap] = calcAutoCorMap(obj, rateMap)
           % This function receives the rateMap as an input and calculates
           % the autocorrelation map of the data using the xcorr2 matlab
           % function. It also normalizes the values
           
           autoCorrelationMap = xcorr2(rateMap);
           v = max(abs(autoCorrelationMap(:)));
           autoCorrelationMap = autoCorrelationMap/v;
            
        end
        
        function [x,y] = findLocationOfMaxValue(obj, map)
           % This functions returns the x,y coordinates of the maximum in
           % the given map (an array)
            
            [v,i] = max(abs(map(:)));
            [x,y] = ind2sub(size(map),i);            
        end
        
    end
    
    
    methods (Access=private)
        
    end
    
end

