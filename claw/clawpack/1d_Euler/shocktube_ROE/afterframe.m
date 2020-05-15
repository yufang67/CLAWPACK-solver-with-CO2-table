axis([0 1 0.5 3.5])

[qref_data,tref] = readamrdata(1,Frame,'./qref/');
if (abs(tref - t) > 1e-8)
  error('Reference data time and current time are not compatible');
end;

hold on;
[qref,xref,p] = plotframe1ez(qref_data,mq,'b-');

h = getlegendinfo;
str1 = sprintf('Computed solution (mx = %d)',mx);
str2 = sprintf('Reference solution (mx = %d)',length(xref));
legend([h,p],str1,str2);

hold off;

clear afterframe;
