function [M, G] = mean_from_harmonics(a_vec, j_0, s0, L)

G = [j_0(:); zeros(length(a_vec)-s0,1)]./(L^3*sqrt(4*pi));
M = G.'*a_vec;