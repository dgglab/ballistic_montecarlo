% a fabFeatureGroup is a collection of fabFeature objects. it enables the
% user to group several objects together and move them in unison. 

classdef fabFeatureGroup < handle
    
    properties
        N_Features
        fabFeatures
        
    end
    
    methods
        %constructor
        function obj = fabFeatureGroup()
            obj.fabFeatures={};
            obj.N_Features=0;
        end
        
        %registering features gives the group the ability to set the
        %position of each fabFeature in the group's coordinate system
        function registerFeature(obj,fabFeature_in,theta,delx,dely)
            obj.N_Features=obj.N_Features+1;
            fabFeature_in.rotate_and_offset(theta,delx,dely);
            obj.fabFeatures{obj.N_Features}=fabFeature_in;

        end
        
        %incorporate another fabFeatureGroup into this fabFeatureGroup
        function addFeatureGroup(obj,FeatureGroup_in)
            for i=1:FeatureGroup_in.N_Features
                obj.registerFeature(FeatureGroup_in.fabFeatures{i},0,0,0);
            end
        end
        
        %rotates all contained features
        function rotate_and_offset(obj,theta,delx,dely)

            for i=1:obj.N_Features
                obj.fabFeatures{i}.rotate_and_offset(theta,delx,dely);
            end          
           
        end
        
        %combines overlapping features of the same layer into a single
        %feature
        function joinFeatures(obj)
            for i=1:obj.N_Features
                layerNums(i)=obj.fabFeatures{i}.LayerNum;
            end

            n=0;
            for i=1:6
                r=find(layerNums==(i-1));
                if ~isempty(r)
                    px=[];py=[];
                    for k=1:length(r)
                        [px,py]=polybool('union',obj.fabFeatures{r(k)}.px,...
                           obj.fabFeatures{r(k)}.py,px,py);
                    end
                    %polybool creates single arrays of px and py, with nan
                    %separating features that don't overlap. below we cut
                    %those appart into distinct features.
                    nEntities=sum(isnan(px))+1;
                    if nEntities>1
                        cutPoints=find(isnan(px));
                    else 
                        cutPoints=[];
                    end
                    cutPoints=[0,cutPoints,length(px)+1];
                            
                    for j=1:nEntities
                        r=(cutPoints(j)+1):(cutPoints(j+1)-1);
                        combinedFeatures{n+1}=fabFeature(px(r),py(r),i-1);
                        n=n+1;
                    end
                end
            end
            
            for i=1:obj.N_Features
                obj.fabFeatures{i}.delete;
            end
            obj.N_Features=0;
            obj.fabFeatures={};
            
            for i=1:n
                obj.registerFeature(combinedFeatures{i},0,0,0);
            end
            
        end

        %plots the fabFeatureGroup
        function plot(obj)
            
            clf; hold on;
            plotcolors={'b','r','k','c','m','y'};
            for i=1:obj.N_Features
                if obj.fabFeatures{i}.LayerNum<5
                    obj.fabFeatures{i}.plot(...
                        plotcolors{obj.fabFeatures{i}.LayerNum+1});
                else
                    obj.fabFeatures{i}.plot(plotcolors{6});
                end
            end
            hold off;
            axis equal;
        end
         
        %produces a DXF file for import into CAD software. 
        function DXF_out(obj,filename)
            addpath('DXFLib');
            d1=dxf_open(filename);
            for i=1:obj.N_Features
                d1.layer=obj.fabFeatures{i}.LayerNum;
                dxf_polyline(d1,obj.fabFeatures{i}.px',obj.fabFeatures{i}.py',...
                    zeros(size(obj.fabFeatures{i}.px')));
            end
            dxf_close(d1);
            
        end
        
                % Grows or shrinks a given polygon
        function grow(obj,pad)
            % Assumes coordinates are given in ccw direction
            for i = 1:obj.N_Features
                obj.fabFeatures{i}.grow(pad);
                
            end

        end
        
        end
            
end
