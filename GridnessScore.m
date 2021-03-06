classdef GridnessScore<handle
    %GRIDNESS With this class one can analyze the gridcells
    %   TODO
    
    
    methods
        
        
        function data = readData(~,path)
            % This function receives a path to a binary file that stores
            % data in the format timestamp,x,y (float32) . It returns a
            % data array with x and y spiking coordinates. It is assumed
            % that the maximum/minimum values in either x or y direction is
            % 1 respectivly -1.
            
            fid = fopen(path,'r');
            tmp = fread(fid,'*float32');
            fclose(fid);
            spikingData = [tmp(1:3:end) tmp(2:3:end) tmp(3:3:end)];
            data = spikingData(:,2:3);
        end
        
        
        function [rateMap] = calcRateMap(obj, data, resolution, filter,sigma)
            % This functions calculates a rate map in the given resolution using
            % the spikingCoordinates data. The size of the map is given by
            % the parameters minxy, maxxy.
            
            % use only positive coordinates, scale and round them so that
            % they are  indexable
            maxi = abs(ceil(max(data(:))));
            mini = abs(floor(min(data(:))));
            if maxi > mini
                mini = maxi;
            end
            data = data + mini;
            data = data / (2*mini);
            data = data * resolution;
            data = round(data);
            
           
            % initialize rateMap
            rateMap = zeros(resolution,resolution);
            
            % fill rateMap according to the times the cell fires in that
            % field
            for i=1:size(data)
                j = data(i,1);
                k = data(i,2);
                if j == 0
                    j = 1;
                end
                if k == 0;
                    k = 1;
                end
                
                if isnan(j) || isnan(k)
                    continue;
                end
                
                rateMap(j,k) = rateMap(j,k) + 1;
            end
            
            
            if filter
                rateMap = imgaussfilt(rateMap,sigma);
            end
            rateMap = rateMap / max(rateMap(:));
 
        end
        
        function [autoCorrelationMap] = calcAutoCorMap(obj, rateMap)
            % This function receives the rateMap as an input and calculates
            % the autocorrelation map of the data using the xcorr2 matlab
            % function. It also normalizes the values
            
            %%% OLD VERSION %%%
            %            autoCorrelationMap = xcorr2(rateMap);
            %            v = max(abs(autoCorrelationMap(:)));
            %            autoCorrelationMap = autoCorrelationMap/v;
            %%%%%%%%%%%%%%%%%%%
            
            n = size(rateMap);
            n = n(1)*n(2);
            rateMap = rateMap - mean(rateMap(:));
            autoCorrelationMap = xcorr2(rateMap) / (n * power(std(rateMap(:)),2));
            
        end
        
        function [x,y] = findLocationOfMaxValue(~, map)
            % This functions returns the x,y coordinates of the maximum in
            % the given map (an array)
            
            [v,i] = max(abs(map(:)));
            [x,y] = ind2sub(size(map),i);
        end
        
        function [gridnessScores] = calcGridnessScores(obj,autoCorrMap,inRad,outRad)
            % This function calculates the gridnessScore at every angle in
            % range [0,180]. The paramters inRad & outRad define the radius
            % of the circle the 6 outer firing fields are assumed to be.
            
            dim = size(autoCorrMap);
            
            
            % look for correct map center
            cntr = obj.findLocationOfMaxValue(autoCorrMap);
            
            ringFilter = zeros(dim(1),dim(1));
            
            for i=1:dim(1)
                for j=1:dim(1)
                    % Check whether point is in radius
                    cntrI = power(cntr - i,2);
                    cntrJ = power(cntr - j,2);
                    dist = cntrI + cntrJ;
                    if power(inRad,2) <= dist && dist <= power(outRad,2)
                        ringFilter(i,j) = 1;
                    end
                    
                end
            end
            
            % build ringMap. May be need to check whether the number of
            % pixel in each field is correct?
            ringMap = autoCorrMap.*ringFilter;
            ringMap = ringMap(cntr - outRad - 1 : cntr + outRad + 1, cntr - outRad - 1 : cntr + outRad + 1);
            [nx,ny] = size(ringMap);
            n = nx * ny;
            
            % rotate one map and calculate relation to original map
            corRot = zeros(180,1);
            for idrot=1:180
                rot = imrotate(ringMap, idrot,'crop');
                A = ringMap;
                B = rot;
                
                %                   corRot(idrot) = corr2(A,B);
                
                
                
                % calculate the correlation according to the fromula in
                % paper Sargolini et al.
                
                dotAA = sum(A.*A);
                dotA0 = sum(A.*ones(nx,ny));
                
                
                dotBB = sum(B.*B);
                dotB0 = sum(B.*ones(nx,ny));
                
                dotAB = sum(A.*B);
                
                corRot(idrot) = (n * sum(dotAB) - sum(dotA0) * sum(dotB0)) / (sqrt(n * sum(dotAA) - power(sum(dotA0),2))*sqrt(n * sum(dotBB)) - power(sum(dotB0),2));
                
                
            end
            gridnessScores = corRot;
            
        end
        
        function gridScore = gridScore(~,gridnessScores)
            % This function takes as input the gridnessScores (output of
            % calcGridnessScore) for every
            % angle in range [0,180] and returns the corresponding gridness
            % score
            
            
            % TODO
            gridScore = min(gridnessScores(60),gridnessScores(120)) - max(gridnessScores(30),max(gridnessScores(90),gridnessScores(150)));
            
            
            
        end
        
        
    end
    
    
    methods (Access=private)
        
    end
    
end

