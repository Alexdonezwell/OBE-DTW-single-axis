function plot_OBEDTW_alignment(q,r,t)
% Visualize the OBEDTW alignment
%    q = [0.86528,0.86904,0.87279,0.87998,0.88885,0.89963,0.91478,0.92993,0.77917,0.62293,0.48612,0.40273,0.51641,0.6301,0.58761,0.85152,0.43589,0.35995,0.39428,0.30555,0.25206,-0.05311,0.07884,0.02596,-0.0884,0.0107,0.06625,-0.18069];
%    q = [0.94078,1.00432,0.99084,0.96139,1.0029,0.90722,0.9041,0.75994,0.64847,0.63292,0.6109,0.58712,0.56587,0.55798,0.55008,0.55394,0.57276,0.59198,0.60874,0.58732,0.55692,0.51305,0.48983,0.51276,0.55116,0.58957,0.53861,0.48455,0.62279,0.55342,0.5314,0.60197,0.6264,0.61596,0.62781,0.65114,0.64945,0.61679,0.6079,0.6299,0.61271,0.6147,0.63592,0.64592,0.98293,0.70187,0.63124,0.68883,0.6974,0.72046,0.77319,0.785,0.7814,0.7725,0.7636,0.75499,0.7522,0.7494,0.7466,0.76598];
          q = [1.866,2.6707,4.1891,1.8981,4.4904,4.1132,3.687,3.4234,2.9866,2.6258,2.1452,1.998,1.7587,1.4847,0.93605,0.63229,0.85788,0.62367,0.17445,-0.070321,0.003772,-0.26771,-0.60562,-0.65062,-0.71323,-0.85487,-0.97704,-0.93326,-0.76257,-0.62565,-0.5141,-0.086319,0.0085245,0.02847,0.19873,0.27892,0.2547,0.26158,0.2416,0.037375,0.26284,0.61129,0.44711,0.50984,0.97645,0.72132,0.83581,0.85279,0.89382,-0.012309,-0.93209,-1.0912,-1.7371,-1.7725,-1.6923,-1.5466,-1.4855,-1.521,-1.4068,-1.4032,-1.6468,-2.0107,-2.228,-2.3482,-2.1267,-1.8776,-1.4764,-1.0569,-0.88614,-0.35537,-0.15366,0.3915,0.42486,1.2434,1.1025,1.6119,1.4782,0.97496,0.42012,-0.12512,-0.089402,-0.93623,-0.97046,-1.3595,-1.4328,-1.3255,-1.199,-0.86797,-0.69657,-0.45848,-0.28243,-0.18794,0.39413,0.27633,0.59814,0.31409,0.24012,0.67152,0.46195,0.893,1.0636,0.0070018,0.19736,-0.13469,-0.10168,0.25671,0.28303,0.13994,0.49232,0.37089,0.40615,0.40562,0.24048,-0.24026,-0.063153,-0.59337,-1.1058,-1.3037,-1.4026,-1.8353,-1.8308,-1.5275,-1.1954,-1.1334,-0.71061,-0.42726,-0.07822,0.13392,0.18604,0.78382,0.75093,0.95413,1.2222,0.77845,0.36883,0.33489,-0.17538,-0.20568,-0.29347,-0.4486,-0.28035,-0.28598,-0.42802,-0.41637,-0.4576,-0.55414,-0.63922,-0.7858,-1.0212,-1.125,-1.1026,-1.1964,-1.2739,-1.3528,-1.3899,-1.3955,-1.5039,-1.4635,-1.5624,-1.4725,-1.418,-1.3785,-1.3552,-1.3434,-1.1581,-1.0909,-1.1415,-1.0702,-0.87751,-0.78621,-0.81716,-0.79726,-0.74368,-0.56401,-0.54596,-0.69759,-0.67075,-0.63327,-0.63999,-0.62083,-0.59657,-0.55553,-0.56246,-0.50117,-0.46397,-0.50331,-0.42977,-0.36596,-0.32845,-0.5728,-0.46922,-0.52587,-0.63343,-0.30213,0.033996,-0.004203,0.030962,0.017628,0.035228,0.10256,-0.15831,-0.15489,-0.28826,-0.11814,-0.040541,-0.032273,0.12784,0.30529,0.13785,0.13348,0.25148,0.1362,0.16891,0.0018185,0.022196,0.15748,0.30712,0.4721,0.34125,0.13511,0.21903,0.39865,0.96456,0.72268,0.75567,0.4745];
%          q = [1,2,3,4,5,4,3,2,1,2,3,4,5,6,7,6,5,4,3,2];

%    r = [0.26547,0.21756,0.43252,0.45569,0.44086,0.40492,0.41481,0.43354,0.51606,0.54993,0.52336,0.48799,0.51221,0.49532,0.53214,0.52122,0.48078,0.62456,0.47958,0.69492,0.77177,0.63438,0.70111,0.78674,0.69739,0.61125,0.54654,0.576664];
%    r = [0.6286,0.60787,0.57149,0.53511,0.49873,0.50795,0.52856,0.54917,0.5625,0.51177,0.37816,0.22171,0.30389,0.33973,0.37325,0.41089,0.35094,0.36273,0.42015,0.27338,0.25641,0.31838,0.33885,0.33459,0.32831,0.37343,0.35882,0.32572,0.34922,0.36011,0.34602,0.40663,0.46995,0.40409,0.40141,0.42912,0.38813,0.39576,0.43472,0.4296,0.42656,0.43219,0.43599,0.32112,0.58706,0.66891,0.56907,0.6173,0.56256,0.71762,0.77202,0.75527,0.83337,0.83618,0.58081,0.45781,0.50588,0.60703,0.60425,0.54321];
          r = [0.29896,0.13382,0.017544,-0.79525,0.055677,0.9053,1.4895,2.4406,3.748,3.7032,3.6858,3.4505,3.0827,2.9336,2.0239,1.6737,1.8677,1.4569,1.1089,0.69335,0.27597,0.019908,-0.226,-0.65486,-1.0135,-1.2104,-1.3392,-1.4462,-1.5712,-1.3595,-0.91565,-0.67807,-0.55187,-0.41525,-0.30661,-0.2253,0.1195,0.21372,-0.026684,-0.16554,0.05539,0.18166,0.032672,0.18156,0.52574,0.51944,0.75448,0.70773,0.64356,0.37791,-0.15472,-0.15725,-0.57388,-0.90722,-1.0105,-1.0527,-1.4041,-1.3999,-1.3654,-1.2734,-1.2424,-1.1541,-1.374,-1.6085,-1.8244,-2.2049,-2.3212,-2.3793,-2.2469,-1.9615,-1.7659,-1.0222,-0.50297,0.5189,1.1269,2.3858,2.6977,2.1224,1.5198,0.79993,0.9473,-0.10079,-0.10049,-0.62949,-0.60108,-0.96502,-1.0957,-0.73393,-0.84277,-0.61788,-0.43946,-0.48666,0.20378,0.074416,0.52287,0.31666,0.29883,0.54684,0.7245,0.52255,0.89146,0.52292,0.4372,0.049101,0.13177,0.40225,0.28279,-0.015875,0.1908,-0.045202,-0.18821,-0.63671,-1.0146,-1.3432,-1.5929,-1.6747,-1.5888,-1.338,-1.2276,-1.0596,-0.37152,-0.32231,0.07854,0.050241,0.056331,-0.28744,-0.13956,0.24317,0.33345,0.56051,0.21831,-0.013123,0.08957,-0.06246,-0.54971,-0.67113,-0.90049,-1.0292,-1.1284,-1.1834,-1.134,-1.1674,-1.1595,-1.1025,-1.0934,-1.0015,-1.0549,-1.1443,-1.2576,-1.2906,-1.1871,-1.2189,-1.2536,-1.2964,-1.3818,-1.3361,-1.3545,-1.3458,-1.3522,-1.2438,-1.1334,-1.0936,-1.0699,-1.0566,-0.99445,-0.88043,-0.93591,-0.88777,-0.87249,-0.87045,-0.83848,-0.78414,-0.7552,-0.70083,-0.70281,-0.70146,-0.60209,-0.56692,-0.54744,-0.44409,-0.52152,-0.38142,-0.50048,-0.43548,-0.51709,-0.55327,-0.60057,-0.55101,-0.6212,-0.62582,-0.53568,-0.33979,-0.31738,-0.36918,-0.20928,-0.11825,0.056292,0.053555,0.071481,0.12347,0.30897,0.33266,0.48413,0.42475,0.3123,0.1731,0.22084,0.19049,0.15658,0.055869,0.23223,0.31825,0.37158,0.25163,0.52811,0.49147,0.51047,0.63224,0.41392,0.51603,0.61113,0.61842,0.59854,0.54284,0.45833,0.53039];
%          r = [2,1,2,3,4,5,4,3,2,1,2,3,4,5,6,7,6,5,4,3];

% disp(numel(q));
% disp(numel(r));
 t=0.05;
%    q = znorm(q);
%    r = znorm(r);

  
[dist, dtwM, path] =OBE_cDTW(q,r,0.1,t);
% /a=path;
 disp(dtwM)
 disp(path)
%  disp(dtwM)
FigHandle = handle(figure);
set(FigHandle, 'Position', [250 250 700 300]);

plot(q, 'k', 'LineWidth', 0.55); 
hold on;
plot(r,'b' ,'LineWidth', 0.55);

 for i =2 : size(path,1)
     plot([path(i,1),path(i,2)],[q(path(i,1)),r(path(i,2))]);
     
 end