function [ CMAP ] = Colormap2d(Nsteps,Msteps, AminBmin, AminBmax, AmaxBmin, AmaxBmax )
%Author - Slavko Rebec (email: srebec@stanford.edu)
%Description - This function creates a custom 2D color map for two ranges of intesity 
%              A and B with resolution Nsteps x Msteps, using provided rgb colors at
%              extremes of A and B

%Inputs
%Asteps - number of color steps for A, higher means finer color resolution
%Bsteps - number of color steps for B, higher means finer color resolution
%AminBmin - 3x1 rgb color for when A is minimal and B is minimal, this is
%           the background color
%AminBmin - 3x1 rgb color for when A is maximal and B is minimal, this is
%           the primary color for A
%AminBmin - 3x1 rgb color for when A is minimal and B is maximal, this is
%           the primary color for B
%AminBmin - 3x1 rgb color for when A is maximal and B is maximal, this is
%           the maximum intensity value

%Output
%CMAP - Color map of size Nsteps x Msteps x 3

rgb = [AminBmin; AminBmax; AmaxBmin; AmaxBmax]';
RGB_IM = cat(3, reshape(rgb(1,:),2,2)',reshape(rgb(2,:),2,2)',reshape(rgb(3,:),2,2)');

[X,Y,Z] = meshgrid(linspace(1,2,Nsteps), linspace(1,2,Msteps), 1:3);
CMAP = interp3(RGB_IM, X, Y,Z);

