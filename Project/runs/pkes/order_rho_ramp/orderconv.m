% order of convergence of implicit euler and generalized runge kutta of
% order 4

clear
close all

% load in all data sets
load ./ref/timeref
load ./ref/powerref
load ./rk4/time1/timerk1
load ./rk4/time1/powerrk1
load ./rk4/time2/timerk2
load ./rk4/time2/powerrk2
load ./rk4/time3/timerk3
load ./rk4/time3/powerrk3
load ./rk4/time4/timerk4
load ./rk4/time4/powerrk4
load ./rk4/time5/timerk5
load ./rk4/time5/powerrk5
load ./rk4/time6/timerk6
load ./rk4/time6/powerrk6
load ./rk4/time7/timerk7
load ./rk4/time7/powerrk7
load ./rk4/time8/timerk8
load ./rk4/time8/powerrk8
% load ./ie1/time02/timeie02
% load ./ie1/time02/powerie02
% load ./ie1/time01/timeie01
% load ./ie1/time01/powerie01
% load ./ie1/time1/timeie1
% load ./ie1/time1/powerie1
% load ./ie1/time2/timeie2
% load ./ie1/time2/powerie2
% load ./ie1/time3/timeie3
% load ./ie1/time3/powerie3
% load ./ie1/time4/timeie4
% load ./ie1/time4/powerie4

% point of interest
tpoi = 10.0;

% get poi indices
idxref = find(abs(timeref - tpoi) < 1.e-8);
idxrk1 = find(abs(timerk1 - tpoi) < 1.e-8);
idxrk2 = find(abs(timerk2 - tpoi) < 1.e-8);
idxrk3 = find(abs(timerk3 - tpoi) < 1.e-8);
idxrk4 = find(abs(timerk4 - tpoi) < 1.e-8);
idxrk5 = find(abs(timerk5 - tpoi) < 1.e-8);
idxrk6 = find(abs(timerk6 - tpoi) < 1.e-8);
idxrk7 = find(abs(timerk7 - tpoi) < 1.e-8);
idxrk8 = find(abs(timerk8 - tpoi) < 1.e-8);
% idxie02 = find(abs(timeie02 - tpoi) < 1.e-8);
% idxie01 = find(abs(timeie01 - tpoi) < 1.e-8);
% idxie1 = find(abs(timeie1 - tpoi) < 1.e-8);
% idxie2 = find(abs(timeie2 - tpoi) < 1.e-8);
% idxie3 = find(abs(timeie3 - tpoi) < 1.e-8);
% idxie4 = find(abs(timeie4 - tpoi) < 1.e-8);

% extract powers
pref = powerref(idxref);
prk1 = powerrk1(idxrk1);
prk2 = powerrk2(idxrk2);
prk3 = powerrk3(idxrk3);
prk4 = powerrk4(idxrk4);
prk5 = powerrk5(idxrk5);
prk6 = powerrk6(idxrk6);
prk7 = powerrk7(idxrk7);
prk8 = powerrk8(idxrk8);
% pie02 = powerie02(idxie02);
% pie01 = powerie01(idxie01);
% pie1 = powerie1(idxie1);
% pie2 = powerie2(idxie2);
% pie3 = powerie3(idxie3);
% pie4 = powerie4(idxie4);

% compute errors to reference
errrk1 = abs(prk1 - pref)/pref*100;
errrk2 = abs(prk2 - pref)/pref*100;
errrk3 = abs(prk3 - pref)/pref*100;
errrk4 = abs(prk4 - pref)/pref*100;
errrk5 = abs(prk5 - pref)/pref*100;
errrk6 = abs(prk6 - pref)/pref*100;
errrk7 = abs(prk7 - pref)/pref*100;
errrk8 = abs(prk8 - pref)/pref*100;
% errie02 = abs(pie02 - pref)/pref*100;
% errie01 = abs(pie01 - pref)/pref*100;
% errie1 = abs(pie1 - pref)/pref*100;
% errie2 = abs(pie2 - pref)/pref*100;
% errie3 = abs(pie3 - pref)/pref*100;
% errie4 = abs(pie4 - pref)/pref*100;

% data for plots
x = [2.0,1.0,0.5,0.25,0.125,0.01,0.001,0.0001];
% yie1 = [errie1,errie2,errie3,errie4];
yrk4 = [errrk1,errrk2,errrk3,errrk4,errrk5,errrk6,errrk7,errrk8];

% % ideal convergences
% fie1 = @(xx) yie1(6)*(xx/x(6)).^1;
frk4 = @(xx) yrk4(8)*(xx/x(8)).^4;
% 
% % plots
% loglog(x,yie1,'--.')
% hold on
% loglog(x,fie1(x),'g-')
% loglog(x,yrk4,'r--.')
% loglog(x,frk4(x),'k-')
% grid on
% 
% % write out gnuplot
% fid = fopen('order.gnuplot','w');
% fprintf(fid,'set terminal pdf dashed\n');
% fprintf(fid,'set output ''order.pdf''\n');
% fprintf(fid,'set key bottom right\n');
% fprintf(fid,'set key box linestyle 1 linecolor rgb ''black''\n');
% fprintf(fid,'set key width -3.0\n');
% fprintf(fid,'set log x\n');
% fprintf(fid,'set log y\n');
% fprintf(fid,'set format y ''%%7.1e''\n');
% fprintf(fid,'set mxtics 10\n');
% fprintf(fid,'set grid x y mxtics\n');
% fprintf(fid,'set xlabel ''Time Step [s]''\n');
% fprintf(fid,'set ylabel ''Difference from Reference [%%]''\n');
% fprintf(fid,'set title ''Constant Reactivity, Reference at 1e-5 s''\n');
% fprintf(fid,'plot ''-'' using 1:2 with points pointtype 5 pointsize 1 linecolor rgb ''blue'' title ''1st Order Implicit Euler'', \\\n');
% fprintf(fid,'     ''-'' using 1:2 with points pointtype 7 pointsize 1 linecolor rgb ''red''  title ''4th Order Gen. Runge Kutta'', \\\n');
% fprintf(fid,'     ''-'' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb ''blue'' title ''Exact 1st Order'', \\\n');
% fprintf(fid,'     ''-'' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb ''red'' title ''Exact 4th Order''\n');
% for i = 1:6
%     fprintf(fid,'%d %d\n',x(i),yie1(i));
% end
% fprintf(fid,'e\n');
% for i = 1:6
%     fprintf(fid,'%d %d\n',x(i),yrk4(i));
% end
% fprintf(fid,'e\n');
% for i = 1:6
%     fprintf(fid,'%d %d\n',x(i),fie1(x(i)));
% end
% fprintf(fid,'e\n');
% for i = 1:6
%     fprintf(fid,'%d %d\n',x(i),frk4(x(i)));
% end
% fprintf(fid,'e\n');