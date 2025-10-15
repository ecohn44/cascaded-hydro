clear; clc; close all;

S_diu = load('unit2_DIU.mat');
S_ddu = load('unit2_DDU.mat');

opts = struct('showBands', true, 'showDelta', false, 'titleTag', 'dry');

plotDIUvDDU(S_diu.X, S_diu.std_hat, ...
                   S_ddu.X, S_ddu.std_hat, ...
                   S_ddu.sysparams, './plots', struct('titleTag','dry'), ...
                   S_diu.U_eff, S_ddu.U_eff);

