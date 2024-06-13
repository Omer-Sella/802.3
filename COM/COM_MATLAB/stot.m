function t_params = stot(s_params)
% p 67 R. Mavaddat. (1996). Network scattering parameter. Singapore: World Scientific.
% ISBN 978-981-02-2305-2. http://books.google.com/?id=287g2NkRYxUC&lpg=PA65&dq=T-parameters+&pg=PA67.
[s11, s12, s21, s22] = deal(s_params(1,1,:), s_params(1,2,:), s_params(2,1,:), s_params(2,2,:));
delta = (s11.*s22-s12.*s21);
s21(s21==0)=eps;
t_params = [1./s21, -s22./s21; s11./s21, -delta./s21];
