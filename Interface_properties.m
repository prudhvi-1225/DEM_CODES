%% Defining Properties for Shear Stress-Strain Curve of each spring interface

G=100*10^6;                    %Shear Modulus
H_1=2*10^6;                    %Hardening Value
%H_1 = 0;
tau_y=4*10^6;                %Yeild Shear Stress value
gamma_u=1.5;                     %Max Allowed Plastic Strain 
gamma_p_s=1;                 %Plastic Strain at start of Softening region
%gamma_p_s=0;
tau_s = tau_y + H_1*gamma_p_s;               %Max Shear Stress Point
G_1 = H_1*G/(G+H_1);                         %Slope of Hardening in Shear stress-strain curve 
G_2 = -tau_s/(gamma_u-(gamma_p_s+tau_s/G));  %Slope of Softening in Shear stress-strain curve
H_2 = G*G_2/(G-G_2);                         %Softening Value
gamma_p_s_max = gamma_p_s - tau_s/H_2;

inx = 0:0.01:gamma_u;
iny = [];
for i = inx
    iny = [iny plotinput(i,tau_y,G,tau_s,G_1,G_2,gamma_u,gamma_p_s)];
end
figure(100)
plot(inx,iny,"LineWidth",2);
ylabel("Shear Stress (\tau) (MPa)")
xlabel("Shear Strain (\gamma)")
grid on;
grid minor;


energy = tau_y^2/(2*G) + tau_s*(gamma_u- (gamma_p_s + tau_s/G))/2 + (tau_s + tau_y)*(gamma_p_s + tau_s/G - tau_y/G)/2;

function y = plotinput(x,ty,g,ts,h1,h2,gu,gps)
if x < ty/g
    y = g*x;
elseif x <= (ts/g + gps)
    y = ty + h1*(x- ty/g);
else
    y = (-ts/(gu - (ts/g + gps)))*(x-gu);
end
end


