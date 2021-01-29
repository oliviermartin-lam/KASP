classdef keckTools < handle
    
   methods (Static)
       
       function [ pistonRamp ] = keck_pistonRamp( pistonMax)
           
           nSeg=36;
           
           c1=[25 26 27 28];
           c2=[24 11 12 13 29];
           c3=[23 10 3 4 14 30];
           c4=[22 9 2 5 15 31];
           c5=[21 8  1 6 16 32];
           c6=[20 7 18 17 33];
           c7=[19 36 35 34];
                                 
           pistonRamp=zeros(1,nSeg);
                      
           pistonRamp(c1)=0*pistonMax/6;
           pistonRamp(c2)=1*pistonMax/6;
           pistonRamp(c3)=2*pistonMax/6;
           pistonRamp(c4)=3*pistonMax/6;
           pistonRamp(c5)=4*pistonMax/6;
           pistonRamp(c6)=5*pistonMax/6;
           pistonRamp(c7)=6*pistonMax/6;
           pistonRamp=pistonRamp';                      
       end
         
        function [out,id_l] = idlToMatlabIndex(s,N,indexValid,flag)
            if nargin <4
                flag = 'slopes';
            end
            % defining sorting grid
            gridMatlab = zeros(N);
            for i=1:N^2
                gridMatlab(i) = i;
            end
            gridIDL    = flip(gridMatlab').*indexValid;
            v          = sort(gridIDL(gridIDL~=0));
            for k=1:length(v(:))
                id_l(k)    = find(gridIDL==v(k));
            end
            nV         = sum(indexValid(:)==1);
            
            
            if strcmp(flag,'slopes')
                % Take x and y parts
                sx  = s(1:2:end,:);
                sy  = s(2:2:end,:);
                % re-sort the s vector
                outx = sx(id_l,:);
                outy = sy(id_l,:);
                %concatenate results
                out  = zeros(2*nV,size(s,2));
                out(1:2:end,:) = outx;
                out(2:2:end,:) = outy;
            elseif strcmp(flag,'actu')
                out  = s(id_l,:);
            elseif strcmp(flag,'modes')
                out = s(:,id_l);
            elseif strcmp(flag,'filter')
                out = s(id_l,id_l);
            else
                out = 0;
            end
        end
   end
                
end
